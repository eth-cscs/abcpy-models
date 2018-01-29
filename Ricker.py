from abcpy.probabilisticmodels import ProbabilisticModel, Hyperparameter
import numpy as np

class Ricker(ProbabilisticModel):
    """Ecological model that describes the observed size of animal population over time
        described in [1].

        [1] S. N. Wood. Statistical inference for noisy nonlinear ecological
        dynamic systems. Nature, 466(7310):1102–1104, Aug. 2010.

        Parameters
        ----------
        parameters: list
            Contains the probabilistic models and hyperparameters from which the model derives.
        n_timestep: int, optional
        Number of timesteps. The default value is 100.
    """
    def __init__(self, parameters, n_timestep=100):
        super(Ricker, self).__init__(parameters)

        self.n_timestep = n_timestep
        #Parameter specifying the dimension of the return values of the distribution.
        self.dimension = n_timestep

    def sample_from_distribution(self, k, rng=np.random.RandomState()):
        timeseries_array = [None] * k
        # Initialize local parameters
        parameters = self.get_parameter_values()
        log_r = parameters[0]
        sigma = parameters[1]
        phi = parameters[2]

        for k in range(0, k):
            # Initialize the time-series
            timeseries_obs_size = np.zeros(shape=(self.n_timestep), dtype=np.float)
            timeseries_true_size = np.ones(shape=(self.n_timestep), dtype=np.float)
            for ind in range(1, self.n_timestep - 1):
                timeseries_true_size[ind] = np.exp(log_r + np.log(timeseries_true_size[ind - 1]) - timeseries_true_size[
                    ind - 1] + sigma * rng.normal(0, 1))
                timeseries_obs_size[ind] = rng.poisson(phi * timeseries_true_size[ind])
            timeseries_array[k] = timeseries_obs_size

        # return an array of objects of type Timeseries
        return [True, np.array(timeseries_array)]

    def _check_parameters_before_sampling(self, parameters):
        return True

    def _check_parameters_fixed(self, parameters):
        return True

    def _check_parameters_at_initialization(self, parameters):
        if(len(parameters)!=3):
            raise ValueError('The provided parameters do not have length 3.')

    def pdf(self, x):
        raise NotImplementedError
