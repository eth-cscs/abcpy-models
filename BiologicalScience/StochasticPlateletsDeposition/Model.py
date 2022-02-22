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

import numpy as np
from abcpy.acceptedparametersmanager import AcceptedParametersManager
from abcpy.inferences import InferenceMethod

class DrawFromPrior(InferenceMethod):
    model = None
    rng = None
    n_samples = None
    backend = None

    n_samples_per_param = None  # this needs to be there otherwise it does not instantiate correctly

    def __init__(self, root_models, backend, seed=None, discard_too_large_values=True):
        self.model = root_models
        self.backend = backend
        self.rng = np.random.RandomState(seed)
        self.discard_too_large_values = discard_too_large_values
        # An object managing the bds objects
        self.accepted_parameters_manager = AcceptedParametersManager(self.model)

    def sample(self, n_samples, n_samples_per_param):
        self.n_samples = n_samples
        self.n_samples_per_param = n_samples_per_param
        self.accepted_parameters_manager.broadcast(self.backend, 1)

        # now generate an array of seeds that need to be different one from the other. One way to do it is the
        # following.
        # Moreover, you cannot use int64 as seeds need to be < 2**32 - 1. How to fix this?
        # Note that this is not perfect; you still have small possibility of having some seeds that are equal. Is there
        # a better way? This would likely not change much the performance
        # An idea would be to use rng.choice but that is too
        seed_arr = self.rng.randint(0, np.iinfo(np.uint32).max, size=n_samples, dtype=np.uint32)
        # check how many equal seeds there are and remove them:
        sorted_seed_arr = np.sort(seed_arr)
        indices = sorted_seed_arr[:-1] == sorted_seed_arr[1:]
        # print("Number of equal seeds:", np.sum(indices))
        if np.sum(indices) > 0:
            # the following can be used to remove the equal seeds in case there are some
            sorted_seed_arr[:-1][indices] = sorted_seed_arr[:-1][indices] + 1
        # print("Number of equal seeds after update:", np.sum(sorted_seed_arr[:-1] == sorted_seed_arr[1:]))
        rng_arr = np.array([np.random.RandomState(seed) for seed in sorted_seed_arr])
        rng_pds = self.backend.parallelize(rng_arr)

        parameters_simulations_pds = self.backend.map(self._sample_parameter, rng_pds)
        parameters_simulations = self.backend.collect(parameters_simulations_pds)
        parameters, simulations = [list(t) for t in zip(*parameters_simulations)]

        parameters = np.squeeze(np.array(parameters))
        simulations = np.squeeze(np.array(simulations))
        #parameters = parameters.reshape((parameters.shape[0], parameters.shape[1]))
        #simulations = simulations.reshape((simulations.shape[0], simulations.shape[2], simulations.shape[3],))

        return parameters, simulations

    def sample_in_chunks(self, n_samples, n_samples_per_param, max_chunk_size=10 ** 4):
        """This splits the data generation in chunks. It is useful when generating large datasets with MPI backend,
        which gives an overflow error due to pickling very large objects."""
        parameters_list = []
        simulations_list = []
        samples_to_sample = n_samples
        while samples_to_sample > 0:
            parameters_part, simulations_part = self.sample(min(samples_to_sample, max_chunk_size), n_samples_per_param)
            samples_to_sample -= max_chunk_size
            parameters_list.append(parameters_part)
            simulations_list.append(simulations_part)
        parameters = np.concatenate(parameters_list)
        simulations = np.concatenate(simulations_list)
        return parameters, simulations

    def _sample_parameter(self, rng, npc=None):
        ok_flag = False

        while not ok_flag:
            self.sample_from_prior(rng=rng)
            theta = self.get_parameters(self.model)
            y_sim = self.simulate(self.n_samples_per_param, rng=rng, npc=npc)

            # if there are no potential infinities there (or if we do not check for those).
            # For instance, Lorenz model may give too large values sometimes (quite rarely).
            if np.sum(np.isinf(np.array(y_sim).astype("float32"))) > 0 and self.discard_too_large_values:
                print("y_sim contained too large values for float32; simulating again.")
            else:
                ok_flag = True

        return theta, y_sim

import numpy as np
from abcpy.acceptedparametersmanager import AcceptedParametersManager
from abcpy.inferences import InferenceMethod
from abcpy.output import Journal

class DrawFromPosterior(InferenceMethod):
    model = None
    rng = None
    n_samples = None
    backend = None

    n_samples_per_param = None  # this needs to be there otherwise it does not instantiate correctly

    def __init__(self, root_models, backend, seed=None, discard_too_large_values=True):
        self.model = root_models
        self.backend = backend
        self.rng = np.random.RandomState(seed)
        self.discard_too_large_values = discard_too_large_values
        # An object managing the bds objects
        self.accepted_parameters_manager = AcceptedParametersManager(self.model)
        self.n_samples_per_param = 1

    def sample(self, journal_file):

        journal = Journal.fromFile(journal_file)
        accepted_parameters = journal.get_accepted_parameters(-1)
        accepted_weights = journal.get_weights(-1)
        n_samples = journal.configuration["n_samples"]

        self.accepted_parameters_manager.broadcast(self.backend, 1)
        # Broadcast Accepted parameters and Accepted weights
        self.accepted_parameters_manager.update_broadcast(self.backend, accepted_parameters=accepted_parameters,
                                                          accepted_weights=accepted_weights)

        seed_arr = self.rng.randint(0, np.iinfo(np.uint32).max, size=n_samples, dtype=np.uint32)
        rng_arr = np.array([np.random.RandomState(seed) for seed in seed_arr])
        index_arr = np.arange(0,n_samples,1)
        data_arr = []
        for i in range(len(rng_arr)):
            data_arr.append([rng_arr[i], index_arr[i]])
        data_pds = self.backend.parallelize(data_arr)

        parameters_simulations_pds = self.backend.map(self._sample_parameter, data_pds)
        parameters_simulations = self.backend.collect(parameters_simulations_pds)
        parameters, simulations = [list(t) for t in zip(*parameters_simulations)]

        parameters = np.squeeze(np.array(parameters))
        simulations = np.squeeze(np.array(simulations))

        return parameters, simulations

    def _sample_parameter(self, data, npc=None):
        if isinstance(data, np.ndarray):
            data = data.tolist()
        rng = data[0]
        index = data[1]
        rng.seed(rng.randint(np.iinfo(np.uint32).max, dtype=np.uint32))

        parameter = self.accepted_parameters_manager.accepted_parameters_bds.value()[index]
        print(parameter)
        parameter_list = [x[0] for x in parameter]
        print(parameter_list)
        self.set_parameters(parameter_list)
        param = self.get_parameters()
        print(param)
        y_sim = self.simulate(n_samples_per_param=1)
        #y_sim = self.model[0].forward_simulate(parameter_list,1)
        return parameter, y_sim

