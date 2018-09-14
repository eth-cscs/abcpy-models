import numpy as np
from abcpy.continuousmodels import Uniform
from Model import FFSimulator

# Define Graphical Model
bufferRatio = Uniform([[0], [3]], name='bufferRatio')
bufferAngle = Uniform([[45], [90]], name='bufferAngle')
kW = Uniform([[0], [10]], name='kW')
kS = Uniform([[0], [10]], name='kS')
kD = Uniform([[0], [10]], name='kD')
decay = Uniform([[0], [1]], name='decay')
diffusion = Uniform([[0], [1]], name='diffusion')
ff = FFSimulator([bufferRatio, bufferAngle, kW, kS, kD, decay, diffusion], name = 'ff')

# Example to Generate Data to check it's correct
fftry = FFSimulator([1.0, 50.0, 1.5, 4.0, 1.0, .1, .25], name = 'fftry')
resultfakeobs1 = fftry.forward_simulate([0.5, 50.0, 1.5, 4.0, 1.0, .1, .25], 10)
resultfakeobs2 = fftry.forward_simulate([1.5, 50.0, 1.5, 4.0, 1.0, .1, .25], 10)

# Check the datasets are different
if np.mean(resultfakeobs1[0])==np.mean(resultfakeobs2[0]):
    print('Datasets are not different!')
  
# Define backend
from abcpy.backends import BackendDummy as Backend
backend = Backend()

# Define Statistics
from abcpy.statistics import Identity
statistics_calculator = Identity(degree=1, cross=False)

# Define distance
from Distance import DistanceType1
distance_calculator = DistanceType1(statistics_calculator)

# Check whether the distance works
if distance_calculator.distance(resultfakeobs1, resultfakeobs1)==distance_calculator.distance(resultfakeobs1, resultfakeobs2):
    print('Something may be wrong with the distance!')

# Define kernel and join the defined kernels
from abcpy.perturbationkernel import MultivariateNormalKernel, JointPerturbationKernel
kernelcontinuous = MultivariateNormalKernel([bufferRatio, bufferAngle, kW, kS, kD, decay, diffusion])
kernel = JointPerturbationKernel([kernelcontinuous])

## APMCABC ##
from abcpy.inferences import APMCABC
sampler = APMCABC([ff], [distance_calculator], backend, kernel, seed = 1)
steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file =4, 1000, 1, 0.1, 0.03, 2, 1.0, None
print('APMCABC Inferring')

# We use resultfakeobs1 as our observed dataset
journal_apmcabc = sampler.sample([resultfakeobs1], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file)
print(journal_apmcabc.posterior_mean())
journal_apmcabc.save('apmcabc_fakeobs1.jrnl')
