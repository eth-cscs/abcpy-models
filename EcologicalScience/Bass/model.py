import numpy as np
from scipy.stats import lognorm
import pyNetLogo
from abcpy.probabilisticmodels import ProbabilisticModel, Continuous, InputConnector

class LogNormal(ProbabilisticModel, Continuous):
    def __init__(self, parameters, name='LogNormal'):
        """
        This class implements a probabilistic model following a Lognormal distribution with mean mu and variance sigma.

        Parameters
        ----------
        parameters: list
            Contains the probabilistic models and hyperparameters from which the model derives.
            The list has two entries: from the first entry mean of the underlying normal distribution and from the second entry variance of the underlying normal
            distribution is derived.
            Note that the second value of the list is strictly greater than 0.

        name: string
            The name that should be given to the probabilistic model in the journal file.
        """

        if not isinstance(parameters, list):
            raise TypeError('Input for LogNormal has to be of type list.')
        if len(parameters) < 2:
            raise ValueError('Input for LogNormal has to be of length 2.')

        input_parameters = InputConnector.from_list(parameters)
        super(LogNormal, self).__init__(input_parameters, name)
        self.visited = False

    def _check_input(self, input_values):
        """
        Returns True if the standard deviation is negative.
        """
        if len(input_values) != 2:
            return False

        if input_values[1] <= 0:
            return False
        return True

    def _check_output(self, parameters):
        """
        Checks parameter values that are given as fixed values.
        """
        return True

    def forward_simulate(self, input_values, k, rng=np.random.RandomState(), mpi_comm=None):
        """
        Samples from a normal distribution using the current values for each probabilistic model from which the model derives.

        Parameters
        ----------
        input_values: list
            List of input parameters, in the same order as specified in the InputConnector passed to the init function
        k: integer
            The number of samples that should be drawn.
        rng: Random number generator
            Defines the random number generator to be used. The default value uses a random seed to initialize the generator.

        Returns
        -------
        list: [np.ndarray]
            A list containing the sampled values as np-array.
        """

        mu = input_values[0]
        sigma = input_values[1]
        result = np.array(rng.lognormal(mu, sigma, k))
        return [np.array([x]).reshape(-1, ) for x in result]

    def get_output_dimension(self):
        return 1
        # Why does the following not work here?
        # return self._dimension

    def pdf(self, input_values, x):
        """
        Calculates the probability density function at point x.
        Commonly used to determine whether perturbed parameters are still valid according to the pdf.

        Parameters
        ----------
        input_values: list
            List of input parameters of the from [mu, sigma]
        x: list
            The point at which the pdf should be evaluated.

        Returns
        -------
        Float:
            The evaluated pdf at point x.
        """

        mu = input_values[0]
        sigma = input_values[1]
        pdf = lognorm(scale=np.exp(mu), s=sigma).pdf(x)
        self.calculated_pdf = pdf
        return pdf

class Bass(ProbabilisticModel, Continuous):
    """Generates time dependent ?????? [1].

        [1] Reference Paper ???

        Parameters
        ----------
        parameters: list
            Contains the parameters: [name of the parameters].
        name: string, optional
            Name of the model
        """

    def __init__(self, parameters, name='Bass'):
        if len(parameters) != 5:
            raise RuntimeError('].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 5:
            raise ValueError('Number of parameters of Bass model must be 5.')

        # Check whether input is from correct domain
        H = input_values[0]
        Am = input_values[1]
        AE = input_values[2]
        PM = input_values[3]
        I = input_values[4]

        ### TODO: Check based on their correct domain !!
        # if theta1 <= 0 or theta2 <= 0:
        # why? this does not make much sense, the parameters of the deterministic part could be smaller than 0
        #    return False
        #
        # if sigma_e < 0 or phi < 0 or phi > 1:
        #     return False
        return True

    def _check_output(self, values):
        # if not isinstance(values[0], np.ndarray):
        #     raise ValueError('Output of the normal distribution is always a number.')
        return True

    def get_output_dimension(self):
        return 87

    def get_number_parameters(self):
        return 5

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        # Extract the input parameters
        H = input_values[0]
        Am = input_values[1]
        AE = input_values[2]
        PM = input_values[3]
        I = input_values[4]

        # Do the actual forward simulation
        vector_of_k_samples = self._run_model(H, Am, AE, PM, I, k)
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result

    def _run_model(self, H, Am, AE, PM, I, k):

        netlogo = pyNetLogo.NetLogoLink(gui=False, netlogo_home='/home/rito/netlogo-5.3.1-64', netlogo_version='5')
        #netlogo = pyNetLogo.NetLogoLink(gui=False, netlogo_home='/home/stats/stsqgk/netlogo-5.3.1-64', netlogo_version='5')
        netlogo.load_model('Bass.nlogo')

        # Setup the model
        netlogo.command("setup")

        # Set up the input values of the parameters
        netlogo.command("set H " + str(H))
        netlogo.command("set Am " + str(Am))
        netlogo.command("set AE " + str(AE))
        netlogo.command("set PM " + str(PM))
        netlogo.command("set I " + str(I * 1e+13))

        ## Simulate and save results
        results = []
        for ind in range(k):
            # spin-up the model
            netlogo.command("spin-up")
            # run the model sufficient no. times
            netlogo.command("go-ABC")
            tmp = []
            # the results are collected within NetLogo
            tmp = tmp + list(netlogo.report("ABC_SSB"))
            tmp = tmp + list(netlogo.report("ABC_Rec"))

            for ind in range(30):
                tmp = tmp + list(netlogo.report("ABC_M"+str(ind)))

            for ind in range(30):
                tmp = tmp + list(netlogo.report("ABC_N" + str(ind)))

            results.append(np.array(tmp))
        netlogo.kill_workspace()
        return results

# parameters = [.75,0.0004,0.01,0.02,3e+13]
# model = Bass(parameters)
# y_simulate = model.forward_simulate(parameters,2)
# print(y_simulate)

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

class DrawFromPosteriorFile(InferenceMethod):
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

    def sample(self, file):

        #journal = Journal.fromFile(journal_file)
        accepted_parameters = np.loadtxt(open(file, "rb"), delimiter=",", skiprows=1)
        accepted_parameters[:,-1] = accepted_parameters[:,-1]/1e+13
        n_samples = accepted_parameters.shape[0]

        self.accepted_parameters_manager.broadcast(self.backend, 1)
        # Broadcast Accepted parameters and Accepted weights
        self.accepted_parameters_manager.update_broadcast(self.backend, accepted_parameters=accepted_parameters)

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

        parameter = self.accepted_parameters_manager.accepted_parameters_bds.value()[index,:]
        print(parameter)
        parameter_list = [x for x in parameter]
        print(parameter_list)
        self.set_parameters(parameter_list)
        parameter = self.get_parameters()
        print(parameter)
        y_sim = self.simulate(n_samples_per_param=1)
        #y_sim = self.model[0].forward_simulate(parameter_list,1)
        return parameter, y_sim
