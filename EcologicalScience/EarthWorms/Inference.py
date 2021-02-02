import numpy as np
from abcpy.statistics import Identity
import logging
# logging.basicConfig(level=logging.DEBUG)
from model import EarthWormGunadi2002, LogNormal

# Define Graphical Model
B_0 = LogNormal([[np.log(967)], [0.3536]], name='B_0')
activation_energy = LogNormal([[np.log(0.25)], [0.3536]], name='activation_energy')
energy_tissue = LogNormal([[np.log(3.6)], [0.3536]], name='energy_tissue')
energy_food = LogNormal([[np.log(10.6)], [0.3536]], name='energy_food')
energy_synthesis = LogNormal([[np.log(3.6)], [0.3536]], name='energy_synthesis')
half_saturation_coeff = LogNormal([[np.log(3.5)], [0.3536]], name='half_saturation_coeff')
max_ingestion_rate = LogNormal([[np.log(.7)], [0.3536]], name='max_ingestion_rate')
mass_birth = LogNormal([[np.log(0.011)], [0.3536]], name='mass_birth')
mass_cocoon = LogNormal([[np.log(0.015)], [0.3536]], name='mass_cocoon')
mass_maximum = LogNormal([[np.log(0.5)], [0.3536]], name='mass_maximum')
mass_sexual_maturity = LogNormal([[np.log(0.25)], [0.3536]], name='mass_sexual_maturity')
growth_constant = LogNormal([[np.log(0.177)], [0.3536]], name='growth_constant')
max_reproduction_rate = LogNormal([[np.log(0.182)], [0.3536]], name='max_reproduction_rate')
speed = LogNormal([[np.log(0.004)], [0.3536]], name='speed')

EarthWorm = EarthWormGunadi2002(
    [B_0, activation_energy, energy_tissue, energy_food, energy_synthesis, half_saturation_coeff,
     max_ingestion_rate, mass_birth, mass_cocoon, mass_maximum, mass_sexual_maturity, growth_constant,
     max_reproduction_rate, speed], name='EarthWormG')

# Example to Generate Data to check it's correct
resultfakeobs1 = EarthWorm.forward_simulate(
    [967, 0.25, 3.6, 10.6, 3.6, 3.5, 0.15, 0.011, 0.015, 0.5, 0.25, 0.177, 0.182, 0.004], 2)
# resultfakeobs2 = EarthWorm.forward_simulate([1, 1, 1, 1, 1, 1, 1, 11, 0.015, 0.5, 0.25, 0.177, 0.182, 0.0002], 2)

# # Define backend
from abcpy.backends import BackendMPI as Backend

backend = Backend()
#
# Define Statistics
statistics_calculator = Identity(degree=1, cross=False)
# #print('# Check whether the statistis works')
# print(statistics_calculator.statistics(resultfakeobs1))
# print(statistics_calculator.statistics(resultfakeobs2))
# Define distance
from abcpy.distances import Euclidean

distance_calculator = Euclidean(statistics_calculator)
# print('# Check whether the distance works')
# print(distance_calculator.distance(resultfakeobs1, resultfakeobs1))
# print(distance_calculator.distance(resultfakeobs1, resultfakeobs2))
#
# Define kernel
from abcpy.perturbationkernel import DefaultKernel

kernel = DefaultKernel([B_0, activation_energy, energy_tissue, energy_food, energy_synthesis, half_saturation_coeff,
                        max_ingestion_rate, mass_birth, mass_cocoon, mass_maximum, mass_sexual_maturity,
                        growth_constant,
                        max_reproduction_rate, speed])

## SABC ##
from abcpy.inferences import SABC

sampler = SABC([EarthWorm], [distance_calculator], backend, kernel, seed=1)

steps, epsilon, n_samples, n_samples_per_param, ar_cutoff, full_output, journal_file = 10, 10000, 500, 1, 0.001, 1, None
print('SABC Inferring')
# We use resultfakeobs1 as our observed dataset
journal_sabc = sampler.sample([resultfakeobs1], steps=steps, epsilon=epsilon, n_samples=n_samples,
                              n_samples_per_param=n_samples_per_param, ar_cutoff=ar_cutoff, full_output=full_output,
                              journal_file=journal_file)
print(journal_sabc.posterior_mean())
journal_sabc.save('sabc_earthworm_fakeobs1.jrnl')

from abcpy.output import Journal

jrnl = Journal.fromFile('sabc_earthworm_fakeobs1.jrnl')
print(jrnl.configuration)
fig, ax = jrnl.plot_ESS()
fig.savefig('ess.pdf')
true_param_value = [967, 0.25, 3.6, 10.6, 3.6, 3.5, 0.15, 0.011, 0.015, 0.5, 0.25, 0.177, 0.182, 0.004]
parameters = ['B_0', 'activation_energy', 'energy_tissue', 'energy_food', 'energy_synthesis', 'half_saturation_coeff',
              'max_ingestion_rate', 'mass_birth', 'mass_cocoon', 'mass_maximum', 'mass_sexual_maturity',
              'growth_constant', 'max_reproduction_rate', 'speed']
jrnl.plot_posterior_distr(path_to_save='posterior_1.pdf', parameters_to_show=parameters[:5],
                          true_parameter_values=true_param_value[:5])
jrnl.plot_posterior_distr(path_to_save='posterior_2.pdf', parameters_to_show=parameters[5:10],
                          true_parameter_values=true_param_value[5:10])
jrnl.plot_posterior_distr(path_to_save='posterior_3.pdf', parameters_to_show=parameters[10:],
                          true_parameter_values=true_param_value[10:])
