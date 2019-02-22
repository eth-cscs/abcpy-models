import copy
import logging
import numpy as np
import pickle, random, string, os

import sys, itertools

from abc import ABCMeta, abstractmethod, abstractproperty
from scipy import optimize

from abcpy.acceptedparametersmanager import *
from abcpy.graphtools import GraphTools
from abcpy.jointapprox_lhd import ProductCombination
from abcpy.jointdistances import LinearCombination
from abcpy.output import Journal
from abcpy.perturbationkernel import DefaultKernel
from abcpy.probabilisticmodels import *
from abcpy.utils import cached


class InferenceMethod(GraphTools, metaclass = ABCMeta):
    """
        This abstract base class represents an inference method.

    """

    def __getstate__(self):
        """Cloudpickle is used with the MPIBackend. This function ensures that the backend itself
        is not pickled
        """
        state = self.__dict__.copy()
        del state['backend']
        return state

    @abstractmethod
    def sample(self):
        """To be overwritten by any sub-class:
        Samples from the posterior distribution of the model parameter given the observed
        data observations.
        """
        raise NotImplementedError

    @abstractproperty
    def model(self):
        """To be overwritten by any sub-class: an attribute specifying the model to be used
        """
        raise NotImplementedError

    @abstractproperty
    def rng(self):
        """To be overwritten by any sub-class: an attribute specifying the random number generator to be used
        """
        raise NotImplementedError

    @abstractproperty
    def backend(self):
        """To be overwritten by any sub-class: an attribute specifying the backend to be used."""
        raise NotImplementedError

    @abstractproperty
    def n_samples(self):
        """To be overwritten by any sub-class: an attribute specifying the number of samples to be generated
        """
        raise NotImplementedError

    @abstractproperty
    def n_samples_per_param(self):
        """To be overwritten by any sub-class: an attribute specifying the number of data points in each simulated         data set."""
        raise NotImplementedError


class BaseMethodsWithKernel(metaclass = ABCMeta):
    """
    This abstract base class represents inference methods that have a kernel.
    """

    @abstractproperty
    def kernel(self):
        """To be overwritten by any sub-class: an attribute specifying the transition or perturbation kernel."""
        raise NotImplementedError

    def perturb(self, column_index, epochs = 10, rng=np.random.RandomState()):
        """
        Perturbs all free parameters, given the current weights.
        Commonly used during inference.

        Parameters
        ----------
        column_index: integer
            The index of the column in the accepted_parameters_bds that should be used for perturbation
        epochs: integer
            The number of times perturbation should happen before the algorithm is terminated

        Returns
        -------
        boolean
            Whether it was possible to set new parameter values for all probabilistic models
        """
        current_epoch = 0

        while current_epoch < epochs:

            # Get new parameters of the graph
            new_parameters = self.kernel.update(self.accepted_parameters_manager, column_index, rng=rng)

            self._reset_flags()

            # Order the parameters provided by the kernel in depth-first search order
            correctly_ordered_parameters = self.get_correct_ordering(new_parameters)

            # Try to set new parameters
            accepted, last_index = self.set_parameters(correctly_ordered_parameters, 0)
            if accepted:
                break
            current_epoch+=1

        if current_epoch == 10:
            return [False]

        return [True, correctly_ordered_parameters]


class BaseLikelihood(InferenceMethod, BaseMethodsWithKernel, metaclass = ABCMeta):
    """
    This abstract base class represents inference methods that use the likelihood.
    """
    @abstractproperty
    def likfun(self):
        """To be overwritten by any sub-class: an attribute specifying the likelihood function to be used."""
        raise NotImplementedError


class BaseDiscrepancy(InferenceMethod, BaseMethodsWithKernel, metaclass = ABCMeta):
    """
    This abstract base class represents inference methods using descrepancy.
    """

    @abstractproperty
    def distance(self):
        """To be overwritten by any sub-class: an attribute specifying the distance function."""
        raise NotImplementedError

class SABC(BaseDiscrepancy, InferenceMethod):
    """
    This base class implements a modified version of Simulated Annealing Approximate Bayesian Computation (SABC) of [1] when the prior is non-informative.

    [1] C. Albert, H. R. Kuensch and A. Scheidegger. A Simulated Annealing Approach to
    Approximate Bayes Computations. Statistics and Computing, (2014).

    Parameters
    ----------
    model : list
        A list of the Probabilistic models corresponding to the observed datasets
    distance : abcpy.distances.Distance
        Distance object defining the distance measure used to compare simulated and observed data sets.
    kernel : abcpy.distributions.Distribution
        Distribution object defining the perturbation kernel needed for the sampling.
    backend : abcpy.backends.Backend
        Backend object defining the backend to be used.
    seed : integer, optional
         Optional initial seed for the random number generator. The default value is generated randomly.
    """

    model = None
    distance = None
    kernel = None
    rng = None

    n_samples = None
    n_samples_per_param = None
    epsilon = None

    smooth_distances_bds = None
    all_distances_bds = None

    backend = None

    def __init__(self, root_models, distances, backend, kernel=None, seed=None):
        self.model = root_models
        # We define the joint Linear combination distance using all the distances for each individual models
        self.distance = LinearCombination(root_models, distances)

        if (kernel is None):

            mapping, garbage_index = self._get_mapping()
            models = []
            for mdl, mdl_index in mapping:
                models.append(mdl)
            kernel = DefaultKernel(models)

        self.kernel = kernel
        self.backend = backend
        self.rng = np.random.RandomState(seed)
        self.logger = logging.getLogger(__name__)

        # these are usually big tables, so we broadcast them to have them once
        # per executor instead of once per task
        self.smooth_distances_bds = None
        self.all_distances_bds = None
        self.accepted_parameters_manager = AcceptedParametersManager(self.model)

        self.simulation_counter = 0


    def sample(self, observations, steps, epsilon, n_samples = 10000, n_samples_per_param = 1, beta = 2, delta = 0.2,
               v = 0.3, ar_cutoff = 0.1, resample = None, n_update = None, full_output=0, journal_file = None):
        """Samples from the posterior distribution of the model parameter given the observed
        data observations.

        Parameters
        ----------
        observations : list
            A list, containing lists describing the observed data sets
        steps : integer
            Number of maximum iterations in the sequential algoritm ("generations")
        epsilon : numpy.float
            A proposed value of threshold to start with.
        n_samples : integer, optional
            Number of samples to generate. The default value is 10000.
        n_samples_per_param : integer, optional
            Number of data points in each simulated data set. The default value is 1.
        beta : numpy.float
            Tuning parameter of SABC, default value is 2.
        delta : numpy.float
            Tuning parameter of SABC, default value is 0.2.
        v : numpy.float, optional
            Tuning parameter of SABC, The default value is 0.3.
        ar_cutoff : numpy.float
            Acceptance ratio cutoff, The default value is 0.1.
        resample: int, optional
            Resample after this many acceptance, The default value is None which takes value inside n_samples
        n_update: int, optional
            Number of perturbed parameters at each step, The default value is None which takes value inside n_samples
        full_output: integer, optional
            If full_output==1, intermediate results are included in output journal.
            The default value is 0, meaning the intermediate results are not saved.
        journal_file: str, optional
            Filename of a journal file to read an already saved journal file, from which the first iteration will start.
            The default value is None.

        Returns
        -------
        abcpy.output.Journal
            A journal containing simulation results, metadata and optionally intermediate results.
        """
        global broken_preemptively
        self.sample_from_prior(rng=self.rng)
        self.accepted_parameters_manager.broadcast(self.backend, observations)
        self.epsilon = epsilon
        self.n_samples = n_samples
        self.n_samples_per_param = n_samples_per_param

        if(journal_file is None):
            journal = Journal(full_output)
            journal.configuration["type_model"] = [type(model).__name__ for model in self.model]
            journal.configuration["type_dist_func"] = type(self.distance).__name__
            journal.configuration["type_kernel_func"] = type(self.kernel)
            journal.configuration["n_samples"] = self.n_samples
            journal.configuration["n_samples_per_param"] = self.n_samples_per_param
            journal.configuration["beta"] = beta
            journal.configuration["delta"] = delta
            journal.configuration["v"] = v
            journal.configuration["ar_cutoff"] = ar_cutoff
            journal.configuration["resample"] = resample
            journal.configuration["n_update"] = n_update
            journal.configuration["full_output"] = full_output
        else:
            journal = Journal.fromFile(journal_file)

        accepted_parameters = None
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

        ## Counter whether broken preemptively
        broken_preemptively = False

        for aStep in range(0, steps):
            self.logger.debug("step {}".format(aStep))
            if(aStep==0 and journal_file is not None):
                accepted_parameters=journal.get_accepted_parameters(-1)
                accepted_weights=journal.get_weights(-1)

                #Broadcast Accepted parameters and Accedpted weights
                self.accepted_parameters_manager.update_broadcast(self.backend, accepted_parameters=accepted_parameters, accepted_weights=accepted_weights)

                kernel_parameters = []
                for kernel in self.kernel.kernels:
                    kernel_parameters.append(
                        self.accepted_parameters_manager.get_accepted_parameters_bds_values(kernel.models))

                #Broadcast Accepted Kernel parameters
                self.accepted_parameters_manager.update_kernel_values(self.backend, kernel_parameters=kernel_parameters)

                new_cov_mats = self.kernel.calculate_cov(self.accepted_parameters_manager)
                accepted_cov_mats = []
                for new_cov_mat in new_cov_mats:
                    if not(new_cov_mat.size == 1):
                        accepted_cov_mats.append(beta * new_cov_mat + 0.0001 * np.trace(new_cov_mat) * np.eye(new_cov_mat.shape[0]))
                    else:
                        accepted_cov_mats.append((beta * new_cov_mat + 0.0001 * new_cov_mat).reshape(1,1))

                # Broadcast Accepted Covariance Matrix
                self.accepted_parameters_manager.update_broadcast(self.backend, accepted_cov_mats=accepted_cov_mats)

            # main SABC algorithm
            # print("INFO: Initialization of SABC")
            seed_arr = self.rng.randint(0, np.iinfo(np.uint32).max, size=int(sample_array[aStep]), dtype=np.uint32)
            rng_arr = np.array([np.random.RandomState(seed) for seed in seed_arr])
            index_arr = self.rng.randint(0, self.n_samples, size=int(sample_array[aStep]), dtype=np.uint32)
            data_arr = []
            for i in range(len(rng_arr)):
                data_arr.append([rng_arr[i], index_arr[i]])
            data_pds = self.backend.parallelize(data_arr)

            # 0: update remotely required variables
            self.logger.info("Broadcasting parameters")
            self.epsilon = epsilon

            # 1: Calculate  parameters
            self.logger.info("Initial accepted parameters")
            params_and_dists_pds = self.backend.map(self._accept_parameter, data_pds)
            params_and_dists = self.backend.collect(params_and_dists_pds)
            new_parameters, filenames, index, acceptance, counter = [list(t) for t in zip(*params_and_dists)]

            # Keeping counter of number of simulations
            for count in counter:
                self.simulation_counter+=count

            #new_parameters = np.array(new_parameters)
            index = np.array(index)
            acceptance = np.array(acceptance)

            # Reading all_distances at Initial step
            if aStep == 0:
                index = np.linspace(0, n_samples - 1, n_samples).astype(int).reshape(n_samples, )
                accept = 0

            # Initialize/Update the accepted parameters and their corresponding distances
            if accepted_parameters is None:
                accepted_parameters = new_parameters
            else:
                for ind in range(len(acceptance)):
                    if acceptance[ind] == 1:
                        accepted_parameters[index[ind]] = new_parameters[ind]

            # 1.5: Update the distance and recompute distances from observed data
            self.logger.info("Updating distance")
            distances = self.distance.distances[0].update(filenames, self.accepted_parameters_manager.observations_bds.value(), self.backend)
            self._update_broadcasts(distances=distances)

            # 2: Compute epsilon
            U = self._average_redefined_distance(distances, epsilon * (1 - delta))
            epsilon = self._schedule(U, v)
            #U = np.mean(distances)
            #epsilon = np.percentile(distances, .1 * 100)
            print(epsilon)
            # 4: Show progress and if acceptance rate smaller than a value break the iteration
            if aStep > 0:
                accept = accept + np.sum(acceptance)
                samples_until = samples_until + sample_array[aStep]
                acceptance_rate = accept / samples_until

                msg = ("updates= {:.2f}, epsilon= {}, u.mean={:e}, acceptance rate: {:.2f}"
                        .format(
                            np.sum(sample_array[1:aStep + 1]) / np.sum(sample_array[1:]) * 100,
                            epsilon, U, acceptance_rate
                            )
                        )
                self.logger.debug(msg)
                if acceptance_rate < ar_cutoff:
                    broken_preemptively = True
                    self.logger.debug("Stopping as acceptance rate is lower than cutoff")
                    break

            # 5: Resampling if number of accepted particles greater than resample
            if accept >= resample and U > 1e-100:
                self.logger.info("Weighted resampling")
                weight = np.exp(-distances * delta / U)
                weight = weight / sum(weight)
                index_resampled = self.rng.choice(np.arange(n_samples, dtype=int), n_samples, replace=1, p=weight)
                accepted_parameters = [accepted_parameters[i] for i in index_resampled]
                distances = distances[index_resampled]

                ## Update U and epsilon:
                # # epsilon = epsilon * (1 - delta)
                # U = np.mean(distances)
                # # epsilon = self._schedule(U, v)
                # epsilon = np.percentile(distances, .1 * 100)
                U = self._average_redefined_distance(distances, epsilon * (1 - delta))
                epsilon = self._schedule(U, v)

                ## Print effective sampling size
                print('Resampling: Effective sampling size: ', 1 / sum(pow(weight / sum(weight), 2)))
                accept = 0
                samples_until = 0

                ## Compute and broadcast accepted parameters, accepted kernel parameters and accepted Covariance matrix
                # Broadcast Accepted parameters and add to journal
                self.accepted_parameters_manager.update_broadcast(self.backend, accepted_weights=accepted_weights, accepted_parameters=accepted_parameters)
                # Compute Accepetd Kernel parameters and broadcast them
                kernel_parameters = []
                for kernel in self.kernel.kernels:
                    kernel_parameters.append(
                        self.accepted_parameters_manager.get_accepted_parameters_bds_values(kernel.models))
                self.accepted_parameters_manager.update_kernel_values(self.backend, kernel_parameters=kernel_parameters)
                # Compute Kernel Covariance Matrix and broadcast it
                new_cov_mats = self.kernel.calculate_cov(self.accepted_parameters_manager)
                accepted_cov_mats = []
                for new_cov_mat in new_cov_mats:
                    if not(new_cov_mat.size == 1):
                        accepted_cov_mats.append(beta * new_cov_mat + 0.0001 * np.trace(new_cov_mat) * np.eye(new_cov_mat.shape[0]))
                    else:
                        accepted_cov_mats.append((beta * new_cov_mat + 0.0001 * new_cov_mat).reshape(1,1))

                self.accepted_parameters_manager.update_broadcast(self.backend, accepted_cov_mats=accepted_cov_mats)

                if (full_output == 1 and aStep<= steps-1):
                    ## Saving intermediate configuration to output journal.
                    print('Saving after resampling')
                    journal.add_accepted_parameters(copy.deepcopy(accepted_parameters))
                    journal.add_weights(copy.deepcopy(accepted_weights))
                    journal.add_distances(copy.deepcopy(distances))
                    names_and_parameters = self._get_names_and_parameters()
                    journal.add_user_parameters(names_and_parameters)
                    journal.number_of_simulations.append(self.simulation_counter)
            else:
                ## Compute and broadcast accepted parameters, accepted kernel parameters and accepted Covariance matrix
                # Broadcast Accepted parameters
                self.accepted_parameters_manager.update_broadcast(self.backend, accepted_weights= accepted_weights, accepted_parameters=accepted_parameters)
                # Compute Accepetd Kernel parameters and broadcast them
                kernel_parameters = []
                for kernel in self.kernel.kernels:
                    kernel_parameters.append(
                        self.accepted_parameters_manager.get_accepted_parameters_bds_values(kernel.models))
                self.accepted_parameters_manager.update_kernel_values(self.backend, kernel_parameters=kernel_parameters)
                # Compute Kernel Covariance Matrix and broadcast it
                new_cov_mats = self.kernel.calculate_cov(self.accepted_parameters_manager)
                accepted_cov_mats = []
                for new_cov_mat in new_cov_mats:
                    if not(new_cov_mat.size == 1):
                        accepted_cov_mats.append(beta * new_cov_mat + 0.0001 * np.trace(new_cov_mat) * np.eye(new_cov_mat.shape[0]))
                    else:
                        accepted_cov_mats.append((beta * new_cov_mat + 0.0001 * new_cov_mat).reshape(1,1))

                self.accepted_parameters_manager.update_broadcast(self.backend, accepted_cov_mats=accepted_cov_mats)

                if (full_output == 1 and aStep <= steps-1):
                    ## Saving intermediate configuration to output journal.
                    journal.add_accepted_parameters(copy.deepcopy(accepted_parameters))
                    journal.add_weights(copy.deepcopy(accepted_weights))
                    journal.add_distances(copy.deepcopy(distances))
                    names_and_parameters = self._get_names_and_parameters()
                    journal.add_user_parameters(names_and_parameters)
                    journal.number_of_simulations.append(self.simulation_counter)

        # Add epsilon_arr, number of final steps and final output to the journal
        # print("INFO: Saving final configuration to output journal.")
        if (full_output == 0) or (full_output ==1 and broken_preemptively and aStep<= steps-1):
            journal.add_accepted_parameters(copy.deepcopy(accepted_parameters))
            journal.add_weights(copy.deepcopy(accepted_weights))
            journal.add_distances(copy.deepcopy(distances))
            self.accepted_parameters_manager.update_broadcast(self.backend,accepted_parameters=accepted_parameters,accepted_weights=accepted_weights)
            names_and_parameters = self._get_names_and_parameters()
            journal.add_user_parameters(names_and_parameters)
            journal.number_of_simulations.append(self.simulation_counter)

        journal.configuration["steps"] = aStep + 1
        journal.configuration["epsilon"] = epsilon

        return journal

    def _smoother_distance(self, distance, old_distance):
        """Smooths the distance using the Equation 14 of [1].

        [1] C. Albert, H. R. Kuensch and A. Scheidegger. A Simulated Annealing Approach to
        Approximate Bayes Computations. Statistics and Computing 0960-3174 (2014).

        Parameters
        ----------
        distance: numpy.ndarray
            Current distance between the simulated and observed data
        old_distance: numpy.ndarray
            Last distance between the simulated and observed data

        Returns
        -------
        numpy.ndarray
            Smoothed distance

        """

        smoothed_distance = np.zeros(shape=(len(distance),))

        for ind in range(0, len(distance)):
            if distance[ind] < np.min(old_distance):
                print('aa')
                smoothed_distance[ind] = (distance[ind] / np.min(old_distance)) / len(old_distance)
            else:
                smoothed_distance[ind] = np.mean(np.array(old_distance) < distance[ind])

        return smoothed_distance

    def _average_redefined_distance(self, distance, epsilon):
        """
        Function to calculate the weighted average of the distance
        Parameters
        ----------
        distance: numpy.ndarray
            Distance between simulated and observed data set
        epsilon: float
            threshold

        Returns
        -------
        numpy.ndarray
            Weighted average of the distance
        """
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

    def _update_broadcasts(self, distances):
        def destroy(bc):
            if bc != None:
                bc.unpersist
                # bc.destroy
        if not distances is None:
            self.distances_bds = self.backend.broadcast(distances)

    # define helper functions for map step
    def _accept_parameter(self, data, npc=None):
        """
        Samples a single model parameter and simulate from it until
        accepted with probabilty exp[-rho(x,y)/epsilon].

        Parameters
        ----------
        seed_and_index: list of two integers
            Initial seed for the random number generator and the index of data to operate on

        Returns
        -------
        numpy.ndarray
            accepted parameter
        """
        if(isinstance(data,np.ndarray)):
            data = data.tolist()
        rng=data[0]
        index=data[1]
        rng.seed(rng.randint(np.iinfo(np.uint32).max, dtype=np.uint32))
        acceptance = 0

        counter = 0

        if self.accepted_parameters_manager.accepted_cov_mats_bds == None:

            while acceptance == 0:
                self.sample_from_prior(rng=rng)
                new_theta = self.get_parameters()
                y_sim = self.simulate(self.n_samples_per_param, rng=rng, npc=npc)
                counter+=1
                distance = self.distance.distance(self.accepted_parameters_manager.observations_bds.value(), y_sim)
                acceptance = rng.binomial(1, np.exp(-distance / self.epsilon), 1)
            acceptance = 1
        else:
            ## Select one arbitrary particle:
            index = rng.choice(self.n_samples, size=1)[0]
            ## Sample proposal parameter and calculate new distance:
            theta = self.accepted_parameters_manager.accepted_parameters_bds.value()[index]

            while True:
                perturbation_output = self.perturb(index, rng=rng)
                if perturbation_output[0] and self.pdf_of_prior(self.model, perturbation_output[1]) != 0:
                    new_theta = perturbation_output[1]
                    break
            y_sim = self.simulate(self.n_samples_per_param, rng=rng, npc=npc)
            counter+=1
            distance = self.distance.distance(self.accepted_parameters_manager.observations_bds.value(), y_sim)

            ## Calculate acceptance probability:
            ratio_prior_prob = self.pdf_of_prior(self.model, perturbation_output[1]) / self.pdf_of_prior(self.model,
                self.accepted_parameters_manager.accepted_parameters_bds.value()[index])
            ratio_likelihood_prob = np.exp((self.distances_bds.value()[index] - distance) / self.epsilon)
            acceptance_prob = ratio_prior_prob * ratio_likelihood_prob

            ## If accepted
            if rng.rand(1) < acceptance_prob:
                acceptance = 1
            else:
                distance = np.inf

        filename = os.getcwd() + '/tmp/' + ''.join([random.choice(string.ascii_letters + string.digits) for n in range(32)])
        with open(filename, 'wb') as fp: pickle.dump([new_theta, y_sim], fp)

        return (new_theta, filename, index, acceptance, counter)