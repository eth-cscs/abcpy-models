import unittest

import numpy as np
from abcpy.continuousmodels import ProbabilisticModel, Continuous, InputConnector
from scipy import integrate


class LotkaVolterra(ProbabilisticModel, Continuous):
    """
    LotkaVolterra model; this has 4 parameters and fixed initial conditions x(t=0)=30. y(t=0)=1.

    We integrate the model over time interval [0,T] (with default T=20) and report 10 evenly spaced points in time for
    both x and y (x are preys, y predators).

    If noise==True, we add some log-normal noise over the different observations (with fixed parameters), in order to
    make the model stochastic.
    """

    def __init__(self, parameters, T=20, n_integration_steps=1000, noise=True, name='LotkaVolterra'):

        self.T = T
        self.n_integration_steps = n_integration_steps
        self.noise = noise
        self.X0 = np.array([30, 1])  # capital X is [x,y]
        self.sigma_lognormal = 0.1
        input_parameters = InputConnector.from_list(parameters)
        super(LotkaVolterra, self).__init__(input_parameters, name)

    def forward_simulate(self, input_values, num_forward_simulations, rng=np.random.RandomState()):
        alpha = input_values[0]
        beta = input_values[1]
        gamma = input_values[2]
        delta = input_values[3]

        # integrate the ODE
        dX_dt = self.define_dX_dt(alpha, beta, gamma, delta)  # define the increment to integrate
        X = self.integrate_ODE(dX_dt)

        # take 10 evenly spaced points in time:
        X = X[self.n_integration_steps // 10 - 1::self.n_integration_steps // 10]

        # now duplicate the above for the required num_forward_simulations
        X = np.stack([X] * num_forward_simulations)

        # add lognormal noise (if noise is required); notice that this can give troubles if the output of the
        # integration is 0 or negative (the ODE solution should not be negative, but numerical integration can reach
        # such numbers). For this reason, if the numbers are negative I put those to a very small positive one:
        if self.noise:
            X[X <= 0] = 1e-20
            X = rng.lognormal(np.log(X), sigma=self.sigma_lognormal)

        return [x for x in X]

    @staticmethod
    def define_dX_dt(alpha, beta, gamma, delta, ):
        def dX_dt(X, t=0):
            """ Return the growth rate of populations. """
            return np.array([alpha * X[0] - beta * X[0] * X[1],
                             -gamma * X[1] + delta * X[0] * X[1]])

        return dX_dt

    def integrate_ODE(self, dX_dt):
        t = np.linspace(0, self.T, self.n_integration_steps)  # time
        X = integrate.odeint(dX_dt, self.X0, t, full_output=False)
        return X

    def get_output_dimension(self):
        return 20

    def _check_input(self, input_values):
        """
        """
        if len(input_values) != 4:
            return False

        # the parameters have to be positive
        if np.any(np.array(input_values) < 0):
            return False

        return True

    def _check_output(self, values):
        return True


class Multivariate_g_and_k_Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.alpha = np.exp(-0.125)
        self.beta = 1
        self.gamma = 0.4
        self.delta = 0.5

        self.model = LotkaVolterra([self.alpha, self.beta, self.gamma, self.delta], noise=True,
                                   n_integration_steps=10000)
        self.model_no_noise = LotkaVolterra([self.alpha, self.beta, self.gamma, self.delta], noise=False,
                                            n_integration_steps=10000)
        self.rng = np.random.RandomState(seed=42)

    def test_check_input(self):
        self.assertTrue(not self.model._check_input([- self.alpha, self.beta, self.gamma, self.delta]))
        self.assertTrue(not self.model._check_input([self.alpha, - self.beta, self.gamma, self.delta]))
        self.assertTrue(not self.model._check_input([self.alpha, self.beta, - self.gamma, self.delta]))
        self.assertTrue(not self.model._check_input([self.alpha, self.beta, self.gamma, - self.delta]))

    def test_forward_sim(self):
        out = self.model.forward_simulate([self.alpha, self.beta, self.gamma, self.delta], num_forward_simulations=2,
                                          rng=self.rng)

        self.assertTrue(not np.allclose(out[0], out[1]))  # these should be different as I am using noise
        self.assertAlmostEqual(np.mean(out[0]), 0.8459084101694787)
        self.assertAlmostEqual(np.std(out[0]), 2.1329435914432353)

        out = self.model_no_noise.forward_simulate([self.alpha, self.beta, self.gamma, self.delta],
                                                   num_forward_simulations=2, rng=self.rng)

        self.assertTrue(np.allclose(out[0], out[1]))  # these should now be equal
        self.assertAlmostEqual(np.mean(out[0]), 0.817750307792599)
        self.assertAlmostEqual(np.std(out[0]), 2.1015960549106514)
