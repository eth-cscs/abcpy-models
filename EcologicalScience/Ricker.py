import numpy as np
from abcpy.continuousmodels import ProbabilisticModel, Continuous, InputConnector
import unittest


class Ricker(ProbabilisticModel, Continuous):
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

    def __init__(self, parameters, n_timestep=100, name="Ricker"):
        input_parameters = InputConnector.from_list(parameters)

        self.n_timestep = n_timestep
        # Parameter specifying the dimension of the return values of the distribution.
        super(Ricker, self).__init__(input_parameters, name)

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        timeseries_array = [None] * k
        # Initialize local parameters
        log_r = input_values[0]
        sigma = input_values[1]
        phi = input_values[2]

        for i in range(0, k):
            # Initialize the time-series
            timeseries_obs_size = np.zeros(shape=self.n_timestep, dtype=np.float)
            timeseries_true_size = np.ones(shape=self.n_timestep, dtype=np.float)
            for ind in range(1, self.n_timestep - 1):
                timeseries_true_size[ind] = np.exp(log_r + np.log(timeseries_true_size[ind - 1]) - timeseries_true_size[
                    ind - 1] + sigma * rng.normal(0, 1))
                timeseries_obs_size[ind] = rng.poisson(phi * timeseries_true_size[ind])
            timeseries_array[i] = timeseries_obs_size

        # return an array of objects of type Timeseries
        return timeseries_array

    def _check_input(self, input_values):
        """
        """
        if len(input_values) != 3:
            return False
        # shold also check range for the different parameters
        return True

    def _check_output(self, values):
        return True

    def get_output_dimension(self):
        return self.n_timestep


class RickerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.log_r = 1
        self.sigma = 1.5
        self.phi = 0.4

        self.model = Ricker([self.log_r, self.sigma, self.phi])
        self.rng = np.random.RandomState(seed=42)

    def test_forward_sim(self):
        out = self.model.forward_simulate([self.log_r, self.sigma, self.phi], k=2, rng=self.rng)
        self.assertEqual(len(np.where(out[0])[0]), 13)
        self.assertTrue(np.allclose(np.where(out[0]), np.array([1, 28, 30, 35, 36, 38, 42, 43, 47, 49, 55, 60, 74])))
