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

        # We expect input of type parameters = [noAP, noNAP, SR_x, pAd, pAg, pT, pF, aT, v_z_AP, v_z_NAP]
        if not isinstance(parameters, list):
            raise TypeError('Input of PlateletDeposition model is of type list')

        if len(parameters) != 10:
            raise RuntimeError('Input list must be of length 10, containing [noAP, noNAP, SR_x, pAd, pAg, pT, pF, aT, v_z_AP, v_z_NAP].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)


    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 10:
            raise ValueError('Number of parameters of PlateletDeposition model must be 10.')

        # Check whether input is from correct domain
        noAP = input_values[0]       # num of Activated PLTs      # PARAMETER
        noNAP = input_values[1]      # num of non-Activated PLTs  # PARAMETER
        SR_x = input_values[2]       # Shear rate                 # PARAMETER
        pAd = input_values[3]        # adhesion rate              # PARAMETER
        pAg = input_values[4]        # aggregation rate           # PARAMETER
        pT = input_values[5]         # deposition rate            # PARAMETER
        pF = input_values[6]         # deposition rate of albumin # PARAMETER
        aT = input_values[7]         # attenuation factor         # PARAMETER
        v_z_AP = input_values[8]     # velocity unit              # PARAMETER
        v_z_NAP = input_values[9]    # velocity unit              # PARAMETER

        if not isinstance(noAP, (int, np.int64, np.int32, np.int16)) or noAP <= 0 or noAP >= 20000:
            print('number of Activated PLTs is not of correct type or out of range')
            return False

        if not isinstance(noNAP, (int, np.int64, np.int32, np.int16)) or noNAP <= 0 or noNAP >= 500000:
            print('number of non-Activated PLTs is not of correct type or out of range')
            return False

        if not isinstance(SR_x, (float, np.float64, np.float32, np.float16)) or SR_x <= 0 or SR_x >= 500:
            print('Shear rate is not of correct type or out of range')
            return False

        if not isinstance(pAd, (float, np.float64, np.float32, np.float16)) or pAd <= 5 or pAd >= 150:
            print('adhesion rate is not of correct type or out of range')
            return False

        if not isinstance(pAg, (float, np.float64, np.float32, np.float16)) or pAg <= 5 or pAg >= 150:
            print('aggregation rate is not of correct type or out of range')
            return False

        if not isinstance(pT, (float, np.float64, np.float32, np.float16)) or pT <= 0.1 or pT >= 10.0:
            print('deposition rate is not of correct type or out of range')
            return False

        if not isinstance(pF, (float, np.float64, np.float32, np.float16)) or pF <= 0.1e-3  or pF >= 9.0e-3:
            print('deposition rate of albumin is not of correct type or out of range')
            return False

        if not isinstance(aT, (float, np.float64, np.float32, np.float16)) or aT <= 0  or aT >= 10:
            print('attenution rate is not of correct type or out of range')
            return False

        if not isinstance(v_z_AP, (float, np.float64, np.float32, np.float16)) or v_z_AP < 1.0e-3  or v_z_AP > 9.0e-3:
            print('v_z_AP is not of correct type or out of range')
            return False

        if not isinstance(v_z_NAP, (float, np.float64, np.float32, np.float16)) or v_z_NAP < 1.0e-4  or v_z_NAP > 9.0e-4:
            print('v_z_NAP is not of correct type or out of range')
            return False

        return True


    def _check_output(self, values):
        if not isinstance(values[0], np.ndarray):
            raise ValueError('Output should be a numpy array.')
        return True

    def get_output_dimension(self):
        return 1

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        if self._check_input(input_values):
            # Extract the input parameters
            noAP = input_values[0]       # num of Activated PLTs      # PARAMETER
            noNAP = input_values[1]      # num of non-Activated PLTs  # PARAMETER
            SR_x = input_values[2]       # Shear rate                 # PARAMETER
            pAd = input_values[3]        # adhesion rate              # PARAMETER
            pAg = input_values[4]        # aggregation rate           # PARAMETER
            pT = input_values[5]         # deposition rate            # PARAMETER
            pF = input_values[6]         # deposition rate of albumin # PARAMETER
            aT = input_values[7]         # attenuation factor         # PARAMETER
            v_z_AP = input_values[8]     # velocity unit              # PARAMETER
            v_z_NAP = input_values[9]    # velocity unit              # PARAMETER
        else:
            raise ValueError('Some values are of not correct')

        # Do the actual forward simulation

        vector_of_k_samples = self.platelet_deposition(noAP, noNAP, SR_x, pAd, pAg, pT, pF, aT, v_z_AP, v_z_NAP, k, rng)
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result

    def platelet_deposition(self, noAP, noNAP, SR_x, pAd, pAg, pT, pF, aT, v_z_AP, v_z_NAP, k, rng):
        seed = rng.randint(np.iinfo(np.int32).max)

        nrow, ncol = 5, 5
        mshape = ncol * nrow
        rshape = nrow * ncol * k
        results = np.reshape(model(rshape, k, int(noAP), int(noNAP), SR_x, pAd, pAg, pT, pF, aT, v_z_AP, v_z_NAP, seed), [k, mshape])
        result = [None] * k
        for ind in range(k):
            result[ind] = results[ind]
        return result
