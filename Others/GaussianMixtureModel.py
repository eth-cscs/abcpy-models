import unittest

import numpy as np
from abcpy.continuousmodels import ProbabilisticModel, Continuous, InputConnector


class BivariateGaussianMixtureModel(ProbabilisticModel, Continuous):
    """We implement here a bivariate gaussian mixture model in which the two components have fixed
     covariance matrices; this has parameters p (mixture ratio) and the means mu0 and mu1, which are the
     2-dimensional means of the two components."""

    def __init__(self, parameters, name='BivariateGaussianMixtureModel'):
        input_parameters = InputConnector.from_list(parameters)
        super(BivariateGaussianMixtureModel, self).__init__(input_parameters, name)

        self.cov0 = np.array([[0.5, -0.3], [-0.3, 0.5]])
        self.cov1 = np.array([[0.25, 0], [0, 0.25]])

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):

        p = input_values[0]
        mu0_0 = input_values[1]
        mu0_1 = input_values[2]
        mu1_0 = input_values[3]
        mu1_1 = input_values[4]
        mu0 = np.array([mu0_0, mu0_1])
        mu1 = np.array([mu1_0, mu1_1])

        # draw the bernoulli probabilities
        bernoulli = rng.binomial(n=1, p=p, size=k)

        result = [None] * k
        for i in range(k):
            if bernoulli[i] == 0:
                result[i] = rng.multivariate_normal(mu0, cov=self.cov0)
            else:
                result[i] = rng.multivariate_normal(mu1, cov=self.cov1)
        # could potentially be improved; not sure it is worth.
        return result

    def get_output_dimension(self):
        return 2

    def _check_input(self, input_values):
        """
        """
        if len(input_values) != 5:
            return False
        p = input_values[0]
        mu0_0 = input_values[1]
        mu0_1 = input_values[2]
        mu1_0 = input_values[3]
        mu1_1 = input_values[4]

        if not (0 <= p <= 1):
            return False

        return True

    def _check_output(self, values):
        return True


class BivariateGaussianMixtureModelTests(unittest.TestCase):
    def setUp(self) -> None:
        self.p = 0.4
        self.mu0_0 = 3
        self.mu0_1 = 0.4
        self.mu1_0 = 12.3
        self.mu1_1 = -4

        self.model = BivariateGaussianMixtureModel([self.p, self.mu0_0, self.mu0_1, self.mu1_0, self.mu1_1])
        self.rng = np.random.RandomState(seed=42)

    def test_check_input(self):
        self.assertTrue(not self.model._check_input([self.p, self.mu0_0, self.mu0_1, self.mu1_0]))
        self.assertTrue(not self.model._check_input([- self.p, self.mu0_0, self.mu0_1, self.mu1_0, self.mu1_1]))

    def test_forward_sim(self):
        out = self.model.forward_simulate([self.p, self.mu0_0, self.mu0_1, self.mu1_0, self.mu1_1], k=2, rng=self.rng)
        self.assertTrue(np.allclose(out[0], np.array([3.07199013, 1.29125853])))
