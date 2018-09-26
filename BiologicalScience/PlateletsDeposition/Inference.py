import numpy as np
from abcpy.continuousmodels import Uniform
from Model import PlateletDeposition

# Define Graphical Model
pAd = Uniform([[50], [150]], name='pAD')
pAg = Uniform([[5], [20]], name='pAg')
pT = Uniform([[0.1], [1.5]], name='pT')
pF = Uniform([[0.5e-3], [3e-3]], name='pF')
aT = Uniform([[0], [10]], name='aT')
PD = PlateletDeposition([pAd, pAg, pT, pF, aT], name = 'PD')

# Observed data
data_obs =  [np.array([ [0,0,0,172200,4808],
[20,1689,26.8,155100,1683],
[60,2004,29.9,149400,0],
[120,1968,31.3,140700,0],
[300,1946,36.6,125800,0]])]

data_obs_1 =  [np.array([ [0,0,0,172200,4808],
[20,1689,26.8,155100,1683],
[60,2004,29.9,149400,0],
[120,1968,31.3,140700,0],
[300,1946,36.6,125801,0]])]

# Example to Generate Data to check it's correct
#PDtry = PlateletDeposition([110.0, 14.6, 0.6, 1.7e-3, 6.0], name = 'PD')
#resultfakeobs1 = PDtry.forward_simulate([110.0, 14.6, 0.6, 1.7e-3, 6.0], 1)

# Define backend
from abcpy.backends import BackendDummy as Backend
backend = Backend()

# Define kernel and join the defined kernels
from abcpy.perturbationkernel import DefaultKernel
kernel = DefaultKernel([pAd, pAg, pT, pF, aT])

# Define Statistics
from Statistics import DepositionStatistics
statistics_calculator = DepositionStatistics(degree=1, cross=False)

# Define distance
from Distance import DepositionDistance
distance_calculator = DepositionDistance(statistics_calculator)

print(distance_calculator.distance(data_obs, data_obs_1))

print('Hello')

## APMCABC ##
from abcpy.inferences import APMCABC
sampler = APMCABC([PD], [distance_calculator], backend, kernel, seed = 1)
steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file =4, 8, 1, 0.1, 0.03, 2, 1.0, None
print('APMCABC Inferring')

# We use resultfakeobs1 as our observed dataset
journal_apmcabc = sampler.sample([data_obs], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file)
print(journal_apmcabc.posterior_mean())
journal_apmcabc.save('apmcabc_fakeobs1.jrnl')