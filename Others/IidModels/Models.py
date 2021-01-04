import numpy as np
from abcpy.continuousmodels import ProbabilisticModel, Continuous, InputConnector, Normal


class IidBeta(ProbabilisticModel, Continuous):
    # this is the ABCpy model for an iid beta sample.
    def __init__(self, parameters, iid_size=1, name='Iid_Beta'):

        self.iid_size = iid_size
        input_parameters = InputConnector.from_list(parameters)
        super(IidBeta, self).__init__(input_parameters, name)

    def forward_simulate(self, input_values, num_forward_simulations, rng=np.random.RandomState()):
        alpha = input_values[0]
        beta = input_values[1]
        result = np.array(rng.beta(alpha, beta, (num_forward_simulations, self.iid_size)))
        return [np.array([x]).reshape(-1, ) for x in result]

    def get_output_dimension(self):
        return self.iid_size

    def _check_input(self, input_values):
        """
        Returns False if the standard deviation is negative.
        """
        if len(input_values) != 2:
            return False

        if input_values[1] <= 0 or input_values[0] <= 0:
            return False
        return True

    def _check_output(self, values):
        return all(values > 0) and all(values < 0)


class IidGamma(ProbabilisticModel, Continuous):
    # this is the ABCpy model for an iid gamma sample.
    def __init__(self, parameters, iid_size=1, name='Iid_Gamma'):

        self.iid_size = iid_size
        input_parameters = InputConnector.from_list(parameters)
        super(IidGamma, self).__init__(input_parameters, name)

    def forward_simulate(self, input_values, num_forward_simulations, rng=np.random.RandomState()):
        k = input_values[0]
        theta = input_values[1]
        result = np.array(rng.gamma(k, theta, (num_forward_simulations, self.iid_size)))
        return [np.array([x]).reshape(-1, ) for x in result]

    def get_output_dimension(self):
        return self.iid_size

    def _check_input(self, input_values):
        """
        Returns True if the standard deviation is negative.
        """
        if len(input_values) != 2:
            return False

        if input_values[1] <= 0 or input_values[0] <= 0:
            return False
        return True

    def _check_output(self, values):
        return all(values > 0)


class IidNormal(Normal):
    # this is the ABCpy model for an iid gaussian sample.
    def __init__(self, parameters, iid_size=1, name='Iid_Normal'):
        self.iid_size = iid_size
        super(IidNormal, self).__init__(parameters, name)

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        mu = input_values[0]
        sigma = input_values[1]
        result = np.array(rng.normal(mu, sigma, (k, self.iid_size)))
        return [np.array([x]).reshape(-1, ) for x in result]

    def get_output_dimension(self):
        return self.iid_size
