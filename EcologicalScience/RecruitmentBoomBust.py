import unittest

import numpy as np
from abcpy.continuousmodels import ProbabilisticModel, Continuous, InputConnector
from abcpy.statistics import Statistics
from scipy.stats import skew, kurtosis  # kurtosis computes the excess kurtosis


class RecruitmentBoomBust(ProbabilisticModel, Continuous):
    """Recruitment, Boom and Bust model as described originally in [1]. It is a discrete stochastic temporal
    model used to represent the fluctuation of population size over time.

    [1] Fasiolo, M., Wood, S. N., Hartig, F., & Bravington, M. V. (2018). An extended empirical saddlepoint
    approximation for intractable likelihoods. Electronic Journal of Statistics, 12(1), 1544-1578.

    Parameters
    ----------
    parameters: list
        Contains the probabilistic models and hyperparameters from which the model derives.
    n_timestep: int, optional
        Number of timesteps. The default value is 250.
    burnin: int, optional
        Number of burnin timesteps to discard transient. The default value is 50.
    """

    def __init__(self, parameters, n_timestep=250, burnin=50, name="RecruitmentBoomBust"):
        input_parameters = InputConnector.from_list(parameters)

        self.n_timestep = n_timestep
        self.burnin = burnin
        # Parameter specifying the dimension of the return values of the distribution.
        super(RecruitmentBoomBust, self).__init__(input_parameters, name)

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        timeseries_array = [None] * k
        # Initialize local parameters
        r = input_values[0]
        kappa = input_values[1]
        alpha = input_values[2]
        beta = input_values[3]

        r_plus_1 = 1 + r

        for i in range(0, k):
            # Initialize the time-series
            timeseries = np.zeros(shape=self.n_timestep + self.burnin + 1, dtype=np.int)
            timeseries[0] = kappa  # initial condition is k
            epsilon_t = rng.poisson(beta, self.n_timestep + self.burnin)
            for ind in range(1, self.n_timestep + self.burnin + 1):
                if timeseries[ind - 1] <= kappa:
                    timeseries[ind] = rng.poisson(timeseries[ind - 1] * r_plus_1) + epsilon_t[ind - 1]
                else:
                    timeseries[ind] = rng.binomial(timeseries[ind - 1], alpha) + epsilon_t[ind - 1]

            timeseries_array[i] = timeseries[self.burnin + 1:]

        return timeseries_array

    def _check_input(self, input_values):
        """
        """
        if len(input_values) != 4:
            return False
        r = input_values[0]
        kappa = input_values[1]
        alpha = input_values[2]
        beta = input_values[3]

        if r <= 0 or beta <= 0 or alpha < 0 or alpha > 1 or kappa < 0:
            return False
        return True

    def _check_output(self, values):
        return True

    def get_output_dimension(self):
        return self.n_timestep


class RecruitmentBoomBustStatistics(Statistics):
    """Statistics for the Recruitment, Boom and Bust model, which was originally described in [1]. The statistics we
    implement here are the ones used in [2]

    [1] Fasiolo, M., Wood, S. N., Hartig, F., & Bravington, M. V. (2018). An extended empirical saddlepoint
    approximation for intractable likelihoods. Electronic Journal of Statistics, 12(1), 1544-1578.

    [2] An, Z., Nott, D. J., & Drovandi, C. (2020). Robust Bayesian synthetic likelihood via a semi-parametric approach.
    Statistics and Computing, 30(3), 543-557.
    """

    def __init__(self):
        pass

    def statistics(self, data):
        if isinstance(data, list):
            if np.array(data).shape == (len(data),):
                if len(data) == 1:
                    data = np.array(data).reshape(1, 1)
                data = np.array(data).reshape(len(data), 1)
            else:
                data = np.concatenate(data).reshape(len(data), -1)
        else:
            raise TypeError('Input data should be of type list, but found type {}'.format(type(data)))

        # need to compute the stepwise differences and the stepwise ratios:
        n_observations = data.shape[0]

        results = np.zeros((n_observations, 12))

        for obs_index in range(n_observations):
            differences = data[obs_index][1:] - data[obs_index][0:-1]
            ratios = (data[obs_index][1:] + 1) / (data[obs_index][0:-1] + 1)

            # now compute mean, variance, skewness and kurtosis; notice that, to be consistent with what is done in
            # An et al. 2020 [2], the kurtosis keeps factor 3 and the variance uses 1 ddof.
            results[obs_index] = [
                np.mean(data[obs_index]), np.var(data[obs_index], ddof=1), skew(data[obs_index]),
                kurtosis(data[obs_index], fisher=False),
                np.mean(ratios), np.var(ratios, ddof=1), skew(ratios), kurtosis(ratios, fisher=False),
                np.mean(differences), np.var(differences, ddof=1), skew(differences),
                kurtosis(differences, fisher=False)]
        return results


class RecruitmentBoomBustTests(unittest.TestCase):
    def setUp(self) -> None:
        self.r = 0.4
        self.kappa = 50
        self.alpha = 0.09
        self.beta = 0.05

        self.model = RecruitmentBoomBust([self.r, self.kappa, self.alpha, self.beta], n_timestep=250)
        self.rng = np.random.RandomState(seed=42)

        self.statistics = RecruitmentBoomBustStatistics()

    def test_forward_sim(self):
        out = self.model.forward_simulate([self.r, self.kappa, self.alpha, self.beta], k=2, rng=self.rng)
        self.assertEqual(out[0].shape[0], self.model.get_output_dimension())
        self.assertEqual(np.sum(out[0]), 3239)  # a fixed value

    def test_wrong_params(self):
        self.assertFalse(self.model._check_input([-self.r, self.kappa, self.alpha, self.beta]))
        self.assertFalse(self.model._check_input([-self.r, -self.kappa, self.alpha, self.beta]))
        self.assertFalse(self.model._check_input([-self.r, self.kappa, -self.alpha, self.beta]))
        self.assertFalse(self.model._check_input([-self.r, self.kappa, 1.3, self.beta]))
        self.assertFalse(self.model._check_input([-self.r, self.kappa, self.alpha, -self.beta]))

    def test_summary(self):
        out = self.model.forward_simulate([self.r, self.kappa, self.alpha, self.beta], k=1, rng=self.rng)
        stats = self.statistics.statistics(out)
        self.assertEqual(np.mean(stats), 54.23283777955975)  # a fixed value
