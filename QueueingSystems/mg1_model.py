import numpy as np
from abcpy.continuousmodels import ProbabilisticModel, Continuous, InputConnector, Normal


class MG1Queue(ProbabilisticModel, Continuous):
    """Simulates a M/G/1 queueing system with Uni[theta1, theta2] service times and Exp(theta3) interarrival times.
    It returns the interdeparture time for the first number_steps steps, assuming the queue starts empty."""

    def __init__(self, parameters, number_steps=5, name='Iid_Beta'):

        self.number_steps = number_steps
        input_parameters = InputConnector.from_list(parameters)
        super(MG1Queue, self).__init__(input_parameters, name)

    def forward_simulate(self, input_values, num_forward_simulations, rng=np.random.RandomState()):
        theta1 = input_values[0]
        theta2 = input_values[1]
        theta3 = input_values[2]

        result = [None] * num_forward_simulations
        for i in range(num_forward_simulations):
            result[i] = self.simulate_mg1(theta1, theta2, theta3, rng)
        return result

    def simulate_mg1(self, theta1, theta2, theta3, rng):
        # use here Lindley equation.
        Y = np.zeros(self.number_steps + 1)
        X = np.zeros(self.number_steps + 1)
        A = np.zeros(self.number_steps + 1)
        inter_dep = np.zeros(self.number_steps + 1)

        for i in range(1, self.number_steps + 1):
            A[i] = rng.exponential(scale=1 / theta3)  # scale is 1/rate
            X[i] = rng.uniform(theta1, theta2)
            Y[i] = np.max([0, Y[i - 1] + X[i - 1] - A[i]])
            # compute the inter-departure times:
            inter_dep[i] = A[i] + X[i] + Y[i] - Y[i - 1] - X[i - 1]
        # print(Y)
        return inter_dep[1:]

    def get_output_dimension(self):
        return self.number_steps

    def _check_input(self, input_values):
        """
        """
        if len(input_values) != 3:
            return False

        if input_values[0] < 0 or input_values[1] <= 0 or input_values[2] <= 0 or input_values[0] > input_values[1]:
            return False
        return True

    def _check_output(self, values):
        return True
