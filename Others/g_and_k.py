import unittest

import numpy as np
from abcpy.continuousmodels import ProbabilisticModel, Continuous, InputConnector


class Multivariate_g_and_k(ProbabilisticModel, Continuous):
    """
    Multivariate g_and_k: a multivariate random normal with mean 0 and covariance matrix with 1 elements on the
    diagonal and rho on the upper and lower secondary diagonal is drawn. Then, the g-and-k transformation is
    applied to obtain the multivariate g-and-k.

    When using one single dimension, this is the standard g-and-k distribution.

    See for instance Section 4.5 of [1] for a discussion of the multivariate g-and-k in an ABC setup.

    This implementation was inspired by the R code at https://github.com/dennisprangle/gk, released under
    the GPL-2.0 License. The corresponding R package is discussed in [2], which also contains details on g-and-k
    distribution.

    [1] Jiang, Bai. "Approximate Bayesian computation with Kullback-Leibler divergence as data discrepancy."
    International Conference on Artificial Intelligence and Statistics. PMLR, 2018.
    [2] Prangle, Dennis. "gk: An R Package for the g-and-k and Generalised g-and-h Distributions."
    arXiv preprint arXiv:1706.06889 (2017), url https://arxiv.org/abs/1706.06889.
    """

    def __init__(self, parameters, size=5, name='Iid_Beta'):

        self.size = size
        self.c = 0.8  # fix this
        input_parameters = InputConnector.from_list(parameters)
        super(Multivariate_g_and_k, self).__init__(input_parameters, name)

    def forward_simulate(self, input_values, num_forward_simulations, rng=np.random.RandomState()):
        A = input_values[0]
        B = input_values[1]
        g = input_values[2]
        k = input_values[3]
        rho = input_values[4]

        result = self.draw_multiv_g_and_k(num_forward_simulations, self.size, A, B, g, k, rho, c=self.c, rng=rng)

        return [x for x in result]

    def draw_multiv_g_and_k(self, n, dim, A, B, g, k, rho, c=0.8, rng=np.random.RandomState()):
        """n is the number of samples, dim is the dimensions"""
        # define the covariance matrix first:
        cov = np.eye(dim) + rho * np.eye(dim, k=1) + rho * np.eye(dim, k=-1)

        z = rng.multivariate_normal(mean=np.zeros(dim), cov=cov, size=n)
        return self.z2gk(z, A, B, g, k, c=c)

    @staticmethod
    def z2gk(z, A, B, g, k, c=0.8):
        """Transform a sample from a standard normal (z) to a g-and-k draw. This assumes z, A, B, g, k and c are
        scalars, while z can be an array."""

        z = np.atleast_1d(z)
        z_squared = z ** 2
        if g == 0:
            term1 = 1
        else:
            term1 = (1 + c * np.tanh(g * z * 0.5))
        term2 = np.zeros_like(z)
        zbig = np.isinf(z_squared)
        zsmall = np.logical_not(zbig)
        term2[zbig] = np.sign(z[zbig]) * np.abs(z[zbig]) ** (1 + 2 * k)
        term2[zsmall] = z[zsmall] * (1 + z_squared[zsmall]) ** k

        return A + B * term1 * term2

    def get_output_dimension(self):
        return self.number_steps

    def _check_input(self, input_values):
        """
        It is standard to take B > 0 and fix c = 0.8, and in this case k ≥ 0 guarantees a proper distribution.
        """
        if len(input_values) != 5:
            return False
        A = input_values[0]
        B = input_values[1]
        g = input_values[2]
        k = input_values[3]
        rho = input_values[4]

        if not (-1 <= rho <= 1):
            return False
        if B <= 0 or k < 0:
            return False

        return True

    def _check_output(self, values):
        return True


class Multivariate_g_and_k_Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.A = 1
        self.B = 3
        self.g = 0.4
        self.k = 2.3
        self.rho = 0.4

        self.model = Multivariate_g_and_k([self.A, self.B, self.g, self.k, self.rho])
        self.rng = np.random.RandomState(seed=42)

    def test_check_input(self):
        self.assertTrue(not self.model._check_input([self.A, self.B, self.g, self.k, -2]))
        self.assertTrue(not self.model._check_input([self.A, -self.B, self.g, self.k, self.rho]))
        self.assertTrue(not self.model._check_input([self.A, self.B, self.g, -self.k, -2]))

    def test_forward_sim(self):
        out = self.model.forward_simulate([self.A, self.B, self.g, self.k, self.rho], num_forward_simulations=2,
                                          rng=self.rng)
        self.assertTrue(np.allclose(out[0], np.array([0.72921481, -8.97429356, 0.77753354, 2.90054443, -16.21483016])))
