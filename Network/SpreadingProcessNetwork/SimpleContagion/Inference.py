from abcpy.models import Model
import numpy as np
from abcpy.statistics import Statistics
from abcpy.distributions import Distribution
from scipy.stats import norm
from abcpy.distances import Distance
import networkx as nx
from abcpy.output import Journal
from scipy import optimize


class SABCDiffusion:
    """This base class implements a modified version of Simulated Annealing Approximate Bayesian Computation (SABC) of [1] when the prior is non-informative.

    [1] C. Albert, H. R. Kuensch and A. Scheidegger. A Simulated Annealing Approach to
    Approximate Bayes Computations. Statistics and Computing, (2014).

    Parameters
    ----------
    model : abcpy.models.Model
        Model object that conforms to the Model class.
    distance : abcpy.distances.Distance
        Distance object that conforms to the Distance class.
    kernel : abcpy.distributions.Distribution
        Distribution object defining the perturbation kernel needed for the sampling
    backend : abcpy.backends.Backend
        Backend object that conforms to the Backend class.
    seed : integer, optional
         Optional initial seed for the random number generator. The default value is generated randomly.
    """

    def __init__(self, model, distance, kernel, backend, seed=None):
        self.model = model
        self.distance = distance
        self.kernel = kernel
        self.backend = backend
        self.rng = np.random.RandomState(seed)

    def sample(self, observations, steps, epsilon, n_samples=10000, n_samples_per_param=1, beta=2, delta=0.2, v=0.3,
               ar_cutoff=0.05, resample=None, n_update=None, adaptcov=1, full_output=0):
        """Samples from the posterior distribution of the model parameter given the observed
        data observations.

        Parameters
        ----------
        observations : numpy.ndarray
            Observed data.
        steps : integer
            Number of maximum iterations in the sequential algoritm ("generations")
        epsilon : numpy.float
            An array of proposed values of epsilon to be used at each steps.
        n_samples : integer, optional
            Number of samples to generate. The default value is 10000.
        n_samples_per_param : integer, optional
            Number of data points in each simulated data set. The default value is 1.
        beta : numpy.float
            Tuning parameter of SABC
        delta : numpy.float
            Tuning parameter of SABC
        v : numpy.float, optional
            Tuning parameter of SABC, The default value is 0.3.
        ar_cutoff : numpy.float
            Acceptance ratio cutoff, The default value is 0.5
        resample: int, optional
            Resample after this many acceptance, The default value if n_samples
        n_update: int, optional
            Number of perturbed parameters at each step, The default value if n_samples
        adaptcov : boolean, optional
            Whether we adapt the covariance matrix in iteration stage. The default value TRUE.
        full_output: integer, optional
            If full_output==1, intermediate results are included in output journal.
            The default value is 0, meaning the intermediate results are not saved.

        Returns
        -------
        abcpy.output.Journal
            A journal containing simulation results, metadata and optionally intermediate results.
        """
        journal = Journal(full_output)
        journal.configuration["type_model"] = type(self.model)
        journal.configuration["type_dist_func"] = type(self.distance)
        journal.configuration["type_kernel_func"] = type(self.kernel)
        journal.configuration["n_samples"] = n_samples
        journal.configuration["n_samples_per_param"] = n_samples_per_param
        journal.configuration["beta"] = beta
        journal.configuration["delta"] = delta
        journal.configuration["v"] = v
        journal.configuration["ar_cutoff"] = ar_cutoff
        journal.configuration["resample"] = resample
        journal.configuration["n_update"] = n_update
        journal.configuration["adaptcov"] = adaptcov
        journal.configuration["full_output"] = full_output

        accepted_parameters = np.zeros(shape=(n_samples, len(self.model.get_parameters())))
        distances = np.zeros(shape=(n_samples,))
        smooth_distances = np.zeros(shape=(n_samples,))
        accepted_weights = np.ones(shape=(n_samples, 1))
        all_distances = None
        accepted_cov_mat = None

        if resample == None:
            resample = n_samples
        if n_update == None:
            n_update = n_samples
        sample_array = np.ones(shape=(steps,))
        sample_array[0] = n_samples
        sample_array[1:] = n_update
        ## Acceptance counter to determine the resampling step
        accept = 0
        samples_until = 0

        # Initialize variables that need to be available remotely
        rc = _RemoteContextSABCDiffusion(self.backend, self.model, self.distance, self.kernel, \
                                          observations, n_samples, n_samples_per_param)

        for aStep in range(0, steps):
            # main SABC algorithm
            # print("INFO: Initialization of SABC")
            seed_arr = self.rng.randint(0, np.iinfo(np.uint32).max, size=int(sample_array[aStep]), dtype=np.uint32)
            seed_pds = self.backend.parallelize(seed_arr)

            # 0: update remotely required variables
            # print("INFO: Broadcasting parameters.")
            rc.epsilon = epsilon
            rc._update_broadcasts(self.backend, accepted_parameters, accepted_cov_mat, smooth_distances, all_distances)

            # 1: Calculate  parameters
            # print("INFO: Initial accepted parameter parameters")
            params_and_dists_pds = self.backend.map(rc._accept_parameter, seed_pds)
            params_and_dists = self.backend.collect(params_and_dists_pds)
            new_parameters, new_distances, new_all_parameters, new_all_distances, index, acceptance = [list(t) for t in
                                                                                                       zip(
                                                                                                           *params_and_dists)]
            new_parameters = np.array(new_parameters)
            new_distances = np.array(new_distances)
            new_all_distances = np.concatenate(new_all_distances)
            index = np.array(index)
            acceptance = np.array(acceptance)
            # print('parallel stuff finsihed')
            # Reading all_distances at Initial step
            if aStep == 0:
                index = np.linspace(0, n_samples - 1, n_samples).astype(int).reshape(n_samples, )
                accept = 0
                all_distances = new_all_distances

            # print(index[acceptance == 1])
            # Initialize/Update the accepted parameters and their corresponding distances
            accepted_parameters[index[acceptance == 1], :] = new_parameters[acceptance == 1, :]
            distances[index[acceptance == 1]] = new_distances[acceptance == 1]

            # 2: Smoothing of the distances
            smooth_distances[index[acceptance == 1]] = self._smoother_distance(distances[index[acceptance == 1]],
                                                                               all_distances)

            # 3: Initialize/Update U, epsilon and covariance of perturbation kernel
            if aStep == 0:
                U = self._avergae_redefined_distance(self._smoother_distance(all_distances, all_distances), epsilon)
            else:
                U = np.mean(smooth_distances)
            epsilon = self._schedule(U, v)
            # if accepted_parameters.shape[1] > 1:
            #    accepted_cov_mat = beta*np.cov(np.transpose(accepted_parameters)) + \
            #    0.0001*np.trace(np.cov(np.transpose(accepted_parameters)))*np.eye(accepted_parameters.shape[1])
            # else:
            accepted_cov_mat = (beta * np.var(np.transpose(accepted_parameters[:, 0])) + 0.0001 * (
            np.var(np.transpose(accepted_parameters[:, 0]))) * np.eye(1))[0][0]

            ## 4: Show progress and if acceptance rate smaller than a value break the iteration

            # print("INFO: Saving intermediate configuration to output journal.")
            if full_output == 1:
                journal.add_parameters(accepted_parameters)
                journal.add_weights(accepted_weights)

            if aStep > 0:
                accept = accept + np.sum(acceptance)
                samples_until = samples_until + sample_array[aStep]
                acceptance_rate = accept / samples_until
                print('updates: ', np.sum(sample_array[1:aStep + 1]) / np.sum(sample_array[1:]) * 100, ' epsilon: ',
                      epsilon, \
                      'u.mean: ', U, 'acceptance rate: ', acceptance_rate)
                if acceptance_rate < ar_cutoff:
                    break

            # 5: Resampling if number of accepted particles greater than resample
            if accept >= resample and U > 1e-100:
                ## Weighted resampling:
                weight = np.exp(-smooth_distances * delta / U)
                weight = weight / sum(weight)
                index_resampled = self.rng.choice(np.arange(n_samples), n_samples, replace=1, p=weight)
                accepted_parameters = accepted_parameters[index_resampled, :]
                smooth_distances = smooth_distances[index_resampled]

                ## Update U and epsilon:
                epsilon = epsilon * (1 - delta)
                U = np.mean(smooth_distances)
                epsilon = self._schedule(U, v)

                ## Print effective sampling size
                print('Resampling: Effective sampling size: ', 1 / sum(pow(weight / sum(weight), 2)))
                accept = 0
                samples_until = 0

        # Add epsilon_arr, number of final steps and final output to the journal
        # print("INFO: Saving final configuration to output journal.")
        if full_output == 0:
            journal.add_parameters(accepted_parameters)
            journal.add_weights(accepted_weights)
        journal.configuration["steps"] = aStep + 1
        journal.configuration["epsilon"] = epsilon
        return journal

    def _smoother_distance(self, distance, old_distance):
        """Smooths the distance using the Equation 14 of [1].

        [1] C. Albert, H. R. Kuensch and A. Scheidegger. A Simulated Annealing Approach to
        Approximate Bayes Computations. Statistics and Computing 0960-3174 (2014).
        """

        smoothed_distance = np.zeros(shape=(len(distance),))

        for ind in range(0, len(distance)):
            if distance[ind] < np.min(old_distance):
                smoothed_distance[ind] = (distance[ind] / np.min(old_distance)) / len(old_distance)
            else:
                smoothed_distance[ind] = np.mean(np.array(old_distance) < distance[ind])

        return smoothed_distance

    def _avergae_redefined_distance(self, distance, epsilon):
        if epsilon == 0:
            U = 0
        else:
            U = np.average(distance, weights=np.exp(-distance / epsilon))

        return (U)

    def _schedule(self, rho, v):
        if rho < 1e-100:
            epsilon = 0
        else:
            fun = lambda epsilon: pow(epsilon, 2) + v * pow(epsilon, 3 / 2) - pow(rho, 2)
            epsilon = optimize.fsolve(fun, rho / 2)

        return (epsilon)


class _RemoteContextSABCDiffusion:
    """
    Contains everything that is sent over the network like broadcast vars and map functions
    """

    def __init__(self, backend, model, distance, kernel, observations, n_samples, n_samples_per_param):
        self.model = model
        self.distance = distance
        self.n_samples = n_samples
        self.n_samples_per_param = n_samples_per_param
        self.epsilon = None
        self.kernel = kernel
        # self._smoother_distance = _smoother_distance

        # these are usually big tables, so we broadcast them to have them once
        # per executor instead of once per task
        self.observations_bds = backend.broadcast(observations)
        self.accepted_parameters_bds = None
        self.accepted_cov_mat_bds = None
        self.smooth_distances_bds = None
        self.all_distances_bds = None

    def _update_broadcasts(self, backend, accepted_parameters, accepted_cov_mat, smooth_distances, all_distances):
        def destroy(bc):
            if bc != None:
                bc.unpersist
                # bc.destroy

        if not accepted_parameters is None:
            self.accepted_parameters_bds = backend.broadcast(accepted_parameters)
        if not accepted_cov_mat is None:
            self.accepted_cov_mat_bds = backend.broadcast(accepted_cov_mat)
        if not smooth_distances is None:
            self.smooth_distances_bds = backend.broadcast(smooth_distances)
        if not all_distances is None:
            self.all_distances_bds = backend.broadcast(all_distances)

    def _smoother_distance_remote(self, distance, old_distance):
        """Smooths the distance using the Equation 14 of [1].

        [1] C. Albert, H. R. Kuensch and A. Scheidegger. A Simulated Annealing Approach to
        Approximate Bayes Computations. Statistics and Computing 0960-3174 (2014).
        """

        smoothed_distance = np.zeros(shape=(len(distance),))

        for ind in range(0, len(distance)):
            if distance[ind] < np.min(old_distance):
                smoothed_distance[ind] = (distance[ind] / np.min(old_distance)) / len(old_distance)
            else:
                smoothed_distance[ind] = np.mean(np.array(old_distance) < distance[ind])

        return smoothed_distance

    # define helper functions for map step
    def _accept_parameter(self, seed):
        """
        Samples a single model parameter and simulate from it until
        accepted with probabilty exp[-rho(x,y)/epsilon].

        :type seed: int
        :rtype: np.array
        :return: accepted parameter
        """
        rng = np.random.RandomState(seed)
        self.model.prior.reseed(rng.randint(np.iinfo(np.uint32).max, dtype=np.uint32))
        self.kernel.reseed(rng.randint(np.iinfo(np.uint32).max, dtype=np.uint32))

        all_parameters = []
        all_distances = []
        index = []
        acceptance = 0

        if self.accepted_cov_mat_bds == None:
            while acceptance == 0:
                self.model.sample_from_prior()
                new_theta = self.model.get_parameters()
                all_parameters.append(self.model.get_parameters())
                y_sim = self.model.simulate(self.n_samples_per_param)
                distance = self.distance.distance(self.observations_bds.value(), y_sim)
                all_distances.append(distance)
                acceptance = rng.binomial(1, np.exp(-distance / self.epsilon), 1)
            acceptance = 1
        else:
            ## Select one arbitrary particle:
            index = rng.choice(self.n_samples, size=1)[0]
            ## Sample proposal parameter and calculate new distance:
            theta = self.accepted_parameters_bds.value()[index, :]
            while True:
                self.kernel.set_parameters([theta[0], self.accepted_cov_mat_bds.value(), theta[1]])
                new_theta = self.kernel.sample(1)[0, :]
                theta_is_accepted = self.model.set_parameters(new_theta)
                if theta_is_accepted and self.model.prior.pdf(self.model.get_parameters()) != 0:
                    break
            y_sim = self.model.simulate(self.n_samples_per_param)
            distance = self.distance.distance(self.observations_bds.value(), y_sim)
            smooth_distance = self._smoother_distance_remote([distance], self.all_distances_bds.value())
            ## Calculate acceptance probability:
            ratio_prior_prob = self.model.prior.pdf(new_theta) / self.model.prior.pdf(
                self.accepted_parameters_bds.value()[index, :])
            ratio_likelihood_prob = np.exp((self.smooth_distances_bds.value()[index] - smooth_distance) / self.epsilon)
            acceptance_prob = ratio_prior_prob * ratio_likelihood_prob
            ## If accepted
            if rng.rand(1) < acceptance_prob:
                acceptance = 1
            else:
                distance = np.inf

        return (new_theta, distance, all_parameters, all_distances, index, acceptance)
