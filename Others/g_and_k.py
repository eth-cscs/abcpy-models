import numpy as np
import unittest
from abcpy.continuousmodels import ProbabilisticModel, Continuous, InputConnector
from scipy.optimize import fsolve
from scipy.stats import multivariate_normal


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

    def __init__(self, parameters, size=5, name='Multivariate_g_and_k'):

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

        result = self._draw_multiv_g_and_k(num_forward_simulations, self.size, A, B, g, k, rho, c=self.c, rng=rng)

        return [x for x in result]

    def _create_cov_matrix(self, dim, rho):
        return np.eye(dim) + rho * np.eye(dim, k=1) + rho * np.eye(dim, k=-1)

    def _draw_multiv_g_and_k(self, n, dim, A, B, g, k, rho, c=0.8, rng=np.random.RandomState()):
        """n is the number of samples, dim is the dimensions"""
        # define the covariance matrix first:

        cov = self._create_cov_matrix(dim, rho)

        z = rng.multivariate_normal(mean=np.zeros(dim), cov=cov, size=n)
        return self._z2gk(z, A, B, g, k, c=c)

    @staticmethod
    def _z2gk(z, A, B, g, k, c=0.8):
        """Transform a sample from a standard normal (z) to a g-and-k draw. This assumes z, A, B, g, k and c are
        scalars, while z can be an array. This is basically the Q function (the inverse of cdf) in the standard
        notation"""

        z = np.atleast_1d(z)
        z_squared = z ** 2
        if g == 0:
            term1 = 1
        else:
            term1 = (1 + c * np.tanh(g * z * 0.5))
        term2 = np.zeros_like(z, dtype=float)
        zbig = np.isinf(z_squared)
        zsmall = np.logical_not(zbig)
        term2[zbig] = np.sign(z[zbig]) * np.abs(z[zbig]) ** (1 + 2 * k)
        term2[zsmall] = z[zsmall] * (1 + z_squared[zsmall]) ** k

        return A + B * term1 * term2

    @staticmethod
    def _Q_log_derivative(z, A, B, g, k, c=0.8):
        """Compute the derivative of the Q function (the inverse cumulative function) evaluated at z."""

        z = np.atleast_1d(z)
        z_squared = z ** 2

        # These are used to correct edge cases
        # (likely to be rare so no need for particularly efficient code)
        zbig = np.isinf(z_squared)
        zsmall = np.logical_not(zbig)
        term1 = np.zeros_like(z, dtype=float)
        if k == 0:
            term1 = 0
        else:
            term1[zsmall] = k * np.log(1 + z_squared[zsmall])
            term1[zbig] = 2 * k * np.log(abs(z[zbig]))
        if g == 0:
            term2 = 1
            term4 = 0
        else:
            gz = g * z
            term2 = 1 + c * np.tanh(gz / 2)
            term4 = c * gz / (2 * np.cosh(gz / 2) ** 2)

        term3 = np.zeros_like(z, dtype=float)
        term3[zsmall] = (1 + (2 * k + 1) * z_squared[zsmall]) / (1 + z_squared[zsmall])
        term3[zbig] = 2 * k + 1
        # term4[ is.infinite(z)] = 0

        return np.log(B) + term1 + np.log(term2 * term3 + term4)

    def get_output_dimension(self):
        return self.size

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

    def logpdf(self, x, input_values):
        """x here can be a single (multivariate) observation (1d array) or a set of observations (2d) array, with
         second index denoting the components"""
        if len(x.shape) not in (1, 2):
            raise RuntimeError("Incorrect number of components for x")

        A = input_values[0]
        B = input_values[1]
        g = input_values[2]
        k = input_values[3]
        rho = input_values[4]

        # first: solve the equation x_j = Q(z_j; theta) for z_j, for all components of x; this requires solving
        # numerically. Check how they do in g and k package.
        z = np.zeros_like(x.flatten())
        for j, x_j in enumerate(x.flatten()):
            func = lambda z: self._z2gk(z, A, B, g, k, self.c) - x_j

            z[j] = fsolve(func, x0=np.array([0]), args=())
        z = z.reshape(x.shape)

        # could also do it with one single call -> faster. Not sure if that works however in the same way, as the
        # scaling may be dependent on the dimension

        # func = lambda z: self._z2gk(z, A, B, g, k, self.c) - x
        # z = fsolve(func, x0=np.zeros_like(x), args=())
        # print(np.allclose(z, z2))

        # second: compute the log pdf by using the multivariate normal pdf and the derivative of the quantile function
        # evaluated in z_j
        cov = self._create_cov_matrix(self.size, rho)
        logpdf_mvn = multivariate_normal.logpdf(z, mean=np.zeros(self.size), cov=cov)

        return np.sum(logpdf_mvn) - np.sum(self._Q_log_derivative(z, A, B, g, k, self.c))


class Univariate_g_and_k(Multivariate_g_and_k):
    """Alias the multivariate one for the univariate case """

    def __init__(self, parameters, name='Univariate_g_and_k'):
        rho = 0  # unused
        super(Univariate_g_and_k, self).__init__(parameters + [rho], size=1, name=name)

    def forward_simulate(self, input_values, num_forward_simulations, rng=np.random.RandomState()):
        rho = 0  # this value is needed only for the parent method but unused with size=1
        return super(Univariate_g_and_k, self).forward_simulate(input_values + [rho], num_forward_simulations,
                                                                rng=rng)


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

    def test_logpdf(self):
        out = self.model.forward_simulate([self.A, self.B, self.g, self.k, self.rho], num_forward_simulations=2,
                                          rng=self.rng)
        logpdf1 = self.model.logpdf(out[0], [self.A, self.B, self.g, self.k, self.rho])
        self.assertAlmostEqual(-17.42923883509208, logpdf1, )

        # check now that computing the logpdf with two observations at once works:
        logpdf2 = self.model.logpdf(out[1], [self.A, self.B, self.g, self.k, self.rho])
        logpdf_joint = self.model.logpdf(np.array(out), [self.A, self.B, self.g, self.k, self.rho])
        self.assertAlmostEqual(logpdf_joint, logpdf1 + logpdf2)

    # I've also tested both the logpdf and the sampling routines against the 'gk' R library (in the 1d case only,
    # as that does not implement higher dimensional cases))
