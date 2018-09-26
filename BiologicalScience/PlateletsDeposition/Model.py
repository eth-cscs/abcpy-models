import numpy as np
from model import model
from abcpy.probabilisticmodels import ProbabilisticModel, Continuous, InputConnector, Discrete

class PlateletDeposition(ProbabilisticModel, Continuous):
    """
        Parameters
        ----------
        parameters: list
            Contains the probabilistic models and hyperparameters from which the model derives.

        """
    def __init__(self, parameters, name='PlateletDeposition'):

        self.maxSimulationTime = 60     # time limit for simulation (s) (in case simulation get stuck due to odd parameter's choice)

        # We expect input of type parameters = [pAd, pAg, pT, pF, aT]
        if not isinstance(parameters, list):
            raise TypeError('Input of PlateletDeposition model is of type list')

        if len(parameters) != 5:
            raise RuntimeError('Input list must be of length 5, containing [pAd, pAg, pT, pF, aT].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)


    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 5:
            raise ValueError('Number of parameters of PlateletDeposition model must be 5.')

        # Check whether input is from correct domain
        pAd = input_values[0]       # adhesion rate             # PARAMETER
        pAg = input_values[1]       # aggregation rate          # PARAMETER
        pT = input_values[2]        # deposition rate           # PARAMETER
        pF = input_values[3]        # deposition rate of albumin# parameter
        aT = input_values[4]        # attenuation factor        # parameter

        if not isinstance(pAd, (float, np.float64, np.float32, np.float16)) or pAd <= 50 or pAd >= 150:
            print('adhesion rate is not of correct type or out of range')
            return False

        if not isinstance(pAg, (float, np.float64, np.float32, np.float16)) or pAg <= 5 or pAg >= 20:
            print('aggregation rate is not of correct type or out of range')
            return False

        if not isinstance(pT, (float, np.float64, np.float32, np.float16)) or pT <= 0.1 or pT >= 1.5:
            print('deposition rate is not of correct type or out of range')
            return False

        if not isinstance(pF, (float, np.float64, np.float32, np.float16)) or pF <= 0.5e-3  or pF >= 3e-3:
            print('deposition rate of albumin is not of correct type or out of range')
            return False

        if not isinstance(aT, (float, np.float64, np.float32, np.float16)) or aT <= 0  or aT >= 10:
            print('attenution rate is not of correct type or out of range')
            return False

        return True


    def _check_output(self, values):
        if not isinstance(values[0], np.ndarray):
            raise ValueError('Output should be a numpy array.')
        return True

    def get_output_dimension(self):
        return 1

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        # Extract the input parameters
        pAd = input_values[0]       # adhesion rate             # PARAMETER
        pAg = input_values[1]       # aggregation rate          # PARAMETER
        pT = input_values[2]        # deposition rate           # PARAMETER
        pF = input_values[3]        # deposition rate of albumin# parameter
        aT = input_values[4]        # attenuation factor        # parameter

        # Do the actual forward simulation
        vector_of_k_samples = self.platelet_deposition(pAd, pAg, pT, pF, aT, k, rng)
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result

    def platelet_deposition(self, pAd, pAg, pT, pF, aT, k, rng):
        seed = rng.randint(np.iinfo(np.int32).max)

        nrow, ncol = 5, 5
        mshape = ncol * nrow
        rshape = nrow * ncol * k
        results = np.reshape(model(rshape, k, pAd, pAg, pT, pF, aT, seed), [k, mshape])
        result = [None] * k
        for ind in range(k):
            result[ind] = results[ind]
        return result
