import numpy as np
from abcpy.continuousmodels import Uniform

from Model import StochLorenz95
from Statistics import HakkarainenLorenzStatistics

# Define Graphical Model
theta1 = Uniform([[0], [1]], name='theta1')
theta2 = Uniform([[0], [1]], name='theta2')
sigma_e = 1
phi = 0.4
lorenz = StochLorenz95([theta1, theta2, sigma_e, phi], time_units=4, n_timestep_per_time_unit=30, name='lorenz')

# Example to Generate Data to check it's correct
resultfakeobs1 = lorenz.forward_simulate([.3, .6, sigma_e, phi], 2)
resultfakeobs2 = lorenz.forward_simulate([.4, .6, sigma_e, phi], 1)

# print('# Check the datasets are different')
print(np.mean(resultfakeobs1[0]))
print(np.mean(resultfakeobs1[1]))

# Define backend
from abcpy.backends import BackendDummy as Backend

backend = Backend()
#
# Define Statistics
statistics_calculator = HakkarainenLorenzStatistics(degree=1, cross=False)

print('# Check whether the statistis works')
print(statistics_calculator.statistics(resultfakeobs1))

#
# Define distance
from abcpy.distances import Euclidean

distance_calculator = Euclidean(statistics_calculator)
#
print('# Check whether the distance works')
print(distance_calculator.distance(resultfakeobs1, resultfakeobs1))
print(distance_calculator.distance(resultfakeobs1, resultfakeobs2))

# Define kernel
from abcpy.perturbationkernel import MultivariateNormalKernel, JointPerturbationKernel

# Join the defined kernels
kernel1 = MultivariateNormalKernel([theta1, theta2])
kernel = JointPerturbationKernel([kernel1])
#
## APMCABC ##
from abcpy.inferences import APMCABC

sampler = APMCABC([lorenz], [distance_calculator], backend, kernel, seed=1)
steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file = 4, 100, 1, 0.1, 0.03, 2, 1.0, None
print('APMCABC Inferring')
# We use resultfakeobs1 as our observed dataset
journal_apmcabc = sampler.sample([resultfakeobs1], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff,
                                 covFactor, full_output, journal_file)
print(journal_apmcabc.posterior_mean())
journal_apmcabc.save('apmcabc_lorenz_fakeobs1.jrnl')
