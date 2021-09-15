import numpy as np
import pyNetLogo
from abcpy.probabilisticmodels import ProbabilisticModel, Continuous, InputConnector

from scipy.stats import lognorm

netlogo_home = "/path/to/netlogo"


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


class EarthWormGunadi2002(ProbabilisticModel, Continuous):
    """Generates time dependent mean mass and hatchlings as described in Gunadi et al. 2002 - Full [1].

        [1] Reference Paper ???

        Parameters
        ----------
        parameters: list
            Contains the parameters: [].
        name: string, optional
            Name of the model
        """

    def __init__(self, parameters, name='EarthWormsGunadi2002'):

        if not isinstance(parameters, list):
            raise TypeError('Input of EarthWorm model is of type list')

        if len(parameters) != 14:
            raise RuntimeError('Input list must be of length 14, containing [B_0, activation_energy, energy_tissue, '
                               'energy_food, energy_synthesis, half_saturation_coeff, '
                               'max_ingestion_rate, mass_birth, mass_cocoon, mass_maximum, '
                               'mass_sexual_maturity, growth_constant, max_reproduction_rate, speed].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 14:
            raise ValueError('Number of parameters of EarthWorm model must be 14.')

        # Check whether input is from correct domain
        # Extract the input parameters
        B_0 = input_values[0]
        activation_energy = input_values[1]
        energy_tissue = input_values[2]
        energy_food = input_values[3]
        energy_synthesis = input_values[4]
        half_saturation_coeff = input_values[5]
        max_ingestion_rate = input_values[6]
        mass_birth = input_values[7]
        mass_cocoon = input_values[8]
        mass_maximum = input_values[9]
        mass_sexual_maturity = input_values[10]
        growth_constant = input_values[11]
        max_reproduction_rate = input_values[12]
        speed = input_values[13]

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
        return 27

    def get_number_parameters(self):
        return 14

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        # Extract the input parameters
        B_0 = input_values[0]
        activation_energy = input_values[1]
        energy_tissue = input_values[2]
        energy_food = input_values[3]
        energy_synthesis = input_values[4]
        half_saturation_coeff = input_values[5]
        max_ingestion_rate = input_values[6]
        mass_birth = input_values[7]
        mass_cocoon = input_values[8]
        mass_maximum = input_values[9]
        mass_sexual_maturity = input_values[10]
        growth_constant = input_values[11]
        max_reproduction_rate = input_values[12]
        speed = input_values[13]

        # Do the actual forward simulation
        vector_of_k_samples = self._run_model(B_0, activation_energy, energy_tissue, energy_food, energy_synthesis,
                                              half_saturation_coeff,
                                              max_ingestion_rate, mass_birth, mass_cocoon, mass_maximum,
                                              mass_sexual_maturity, growth_constant,
                                              max_reproduction_rate, speed, k)
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result

    # The following functions define running of NetLogo models
    # Run model 1
    def _run_model(self, B_0, activation_energy, energy_tissue, energy_food, energy_synthesis, half_saturation_coeff,
                   max_ingestion_rate, mass_birth, mass_cocoon, mass_maximum, mass_sexual_maturity, growth_constant,
                   max_reproduction_rate, speed, k):
        ## Initiate netlogo model
        netlogo = pyNetLogo.NetLogoLink(gui=False, netlogo_home=netlogo_home, netlogo_version='5')
        netlogo.load_model('src/models/g2002_1.nlogo')

        ## Set parameter values
        netlogo.command("set B_0 " + str(B_0))
        netlogo.command("set activation_energy " + str(activation_energy))
        netlogo.command("set energy_tissue " + str(energy_tissue))
        netlogo.command("set energy_food " + str(energy_food))
        netlogo.command("set energy_synthesis " + str(energy_synthesis))
        netlogo.command("set half_saturation_coef " + str(half_saturation_coeff))
        netlogo.command("set max_ingestion_rate " + str(max_ingestion_rate))
        netlogo.command("set mass_birth " + str(mass_birth))
        netlogo.command("set mass_cocoon " + str(mass_cocoon))
        netlogo.command("set mass_maximum " + str(mass_maximum))
        netlogo.command("set mass_sexual_maturity " + str(mass_sexual_maturity))
        netlogo.command("set growth_constant " + str(growth_constant))
        netlogo.command("set max_reproduction_rate " + str(max_reproduction_rate))
        netlogo.command("set speed " + str(speed))

        ## Simulate and save results
        result = []
        for ind in range(k):
            result_tmp = np.zeros(shape=(27))
            netlogo.command("setup")
            result_tmp[0] = netlogo.report("mean-mass")
            for t in range(1, 27):
                netlogo.command("repeat 7 [go]")
                result_tmp[t] = netlogo.report("mean-mass")
            result.append(result_tmp.flatten())
        netlogo.kill_workspace()
        return result


class EarthWormGunadiEdwards2003(ProbabilisticModel, Continuous):
    """Generates time dependent mean mass and hatchlings as described in Gunadi & Edwards 2003 - Full [1].

        [1] Reference Paper ???

        Parameters
        ----------
        parameters: list
            Contains the parameters: [].
        name: string, optional
            Name of the model
        """

    def __init__(self, parameters, name='EarthWormsGunadi2002'):

        if not isinstance(parameters, list):
            raise TypeError('Input of EarthWorm model is of type list')

        if len(parameters) != 14:
            raise RuntimeError('Input list must be of length 14, containing [B_0, activation_energy, energy_tissue, '
                               'energy_food, energy_synthesis, half_saturation_coeff, '
                               'max_ingestion_rate, mass_birth, mass_cocoon, mass_maximum, '
                               'mass_sexual_maturity, growth_constant, max_reproduction_rate, speed].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 14:
            raise ValueError('Number of parameters of EarthWorm model must be 14.')

        # Check whether input is from correct domain
        # Extract the input parameters
        B_0 = input_values[0]
        activation_energy = input_values[1]
        energy_tissue = input_values[2]
        energy_food = input_values[3]
        energy_synthesis = input_values[4]
        half_saturation_coeff = input_values[5]
        max_ingestion_rate = input_values[6]
        mass_birth = input_values[7]
        mass_cocoon = input_values[8]
        mass_maximum = input_values[9]
        mass_sexual_maturity = input_values[10]
        growth_constant = input_values[11]
        max_reproduction_rate = input_values[12]
        speed = input_values[13]

        ### TODO: Check based on their correct domain !!
        if B_0 <= 0 or activation_energy < 0 or energy_tissue < 0 or energy_food < 0 or energy_synthesis < 0 or half_saturation_coeff < 0 \
                or max_ingestion_rate < 0 or mass_birth < 0 or mass_cocoon < 0 or mass_maximum < 0 or mass_sexual_maturity < 0 or growth_constant < 0 \
                or max_reproduction_rate < 0 or speed < 0:
            return False

        return True

    def _check_output(self, values):
        # if not isinstance(values[0], np.ndarray):
        #     raise ValueError('Output of the normal distribution is always a number.')
        return True

    def get_output_dimension(self):
        return 122

    def get_number_parameters(self):
        return 14

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        # Extract the input parameters
        B_0 = input_values[0]
        activation_energy = input_values[1]
        energy_tissue = input_values[2]
        energy_food = input_values[3]
        energy_synthesis = input_values[4]
        half_saturation_coeff = input_values[5]
        max_ingestion_rate = input_values[6]
        mass_birth = input_values[7]
        mass_cocoon = input_values[8]
        mass_maximum = input_values[9]
        mass_sexual_maturity = input_values[10]
        growth_constant = input_values[11]
        max_reproduction_rate = input_values[12]
        speed = input_values[13]

        # Do the actual forward simulation
        vector_of_k_samples = self._run_model(B_0, activation_energy, energy_tissue, energy_food, energy_synthesis,
                                              half_saturation_coeff,
                                              max_ingestion_rate, mass_birth, mass_cocoon, mass_maximum,
                                              mass_sexual_maturity, growth_constant,
                                              max_reproduction_rate, speed, k)
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result

    # The following functions define running of NetLogo models
    # Run model 1
    def _run_model(self, B_0, activation_energy, energy_tissue, energy_food, energy_synthesis, half_saturation_coeff,
                   max_ingestion_rate, mass_birth, mass_cocoon, mass_maximum, mass_sexual_maturity, growth_constant,
                   max_reproduction_rate, speed, k):

        ## Initiate netlogo model
        netlogo = pyNetLogo.NetLogoLink(gui=False, netlogo_home=netlogo_home, netlogo_version='5')
        netlogo.load_model('src/models/ge2003_1.nlogo')

        ## Set parameter values
        netlogo.command("set B_0 " + str(B_0))
        netlogo.command("set activation_energy " + str(activation_energy))
        netlogo.command("set energy_tissue " + str(energy_tissue))
        netlogo.command("set energy_food " + str(energy_food))
        netlogo.command("set energy_synthesis " + str(energy_synthesis))
        netlogo.command("set half_saturation_coef " + str(half_saturation_coeff))
        netlogo.command("set max_ingestion_rate " + str(max_ingestion_rate))
        netlogo.command("set mass_birth " + str(mass_birth))
        netlogo.command("set mass_cocoon " + str(mass_cocoon))
        netlogo.command("set mass_maximum " + str(mass_maximum))
        netlogo.command("set mass_sexual_maturity " + str(mass_sexual_maturity))
        netlogo.command("set growth_constant " + str(growth_constant))
        netlogo.command("set max_reproduction_rate " + str(max_reproduction_rate))
        netlogo.command("set speed " + str(speed))

        ## Simulate and save results
        result = []
        for ind in range(k):
            result_tmp = np.zeros(shape=(61, 2))
            netlogo.command("setup")
            result_tmp[0] = netlogo.report("mean-mass")
            for t in range(1, 61):
                netlogo.command("repeat 7 [go]")
                result_tmp[t] = netlogo.report("mean-mass")
            result.append(result_tmp.flatten())
        netlogo.kill_workspace()
        return result

# parameters = [967, 0.25, 7, 10.6, 3.6, 3.5, 0.15, 0.011, 0.015, 0.5, 0.25, 0.177, 0.182, 0.0004]
#
# model = EarthWormGunadi2002(parameters)
# y_simulate = model.forward_simulate(parameters,2)
# print(y_simulate)


# model = EarthWormGunadiEdwards2003(parameters)
# y_simulate = model.forward_simulate(parameters,2)
# print(y_simulate)
