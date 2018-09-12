import numpy as np
from abcpy.continuousmodels import Uniform
from abcpy.discretemodels import Binomial
from Model import FFSimulator, DiscreteUniform

# Some comment
# Define Graphical Model
kW = Uniform([[0], [10]], name='kW')
kS = Uniform([[0], [10]], name='kS')
kD = Uniform([[0], [10]], name='kD')
decay = Uniform([[0], [1]], name='decay')
diffusion = Uniform([[0], [1]], name='diffusion')
bufferRatio = Uniform([[0], [10]], name='bufferratio')
bufferAngle = DiscreteUniform([[45],[90]], name='bufferAngle')
ff = FFSimulator([bufferRatio, bufferAngle, kW, kS, kD, decay, diffusion], name = 'ff')

# Example to Generate Data to check it's correct
fftry = FFSimulator([8.2, 50, 3.2, 4.1, 1.1, .3, .1], name = 'fftry')
resultfakeobs1 = fftry.forward_simulate([8.2, 50, 3.2, 4.1, 1.1, .3, .1], 10)
resultfakeobs2 = fftry.forward_simulate([9.2, 50, 3.2, 4.1, 1.1, .3, .1], 10)

# Check the datasets are different
if np.mean(resultfakeobs1[0])==np.mean(resultfakeobs2[0]):
    print('Datasets are not different!')
  
# Define backend
from abcpy.backends import BackendMPI as Backend
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

# Define kernel
from abcpy.perturbationkernel import MultivariateNormalKernel, RandomWalkKernel, JointPerturbationKernel
# Join the defined kernels
kernelcontinuous = MultivariateNormalKernel([bufferRatio, kW, kS, kD, decay, diffusion])
kerneldiscrete = RandomWalkKernel([bufferAngle])
kernel = JointPerturbationKernel([kernelcontinuous, kerneldiscrete])
#
## APMCABC ##
from abcpy.inferences import APMCABC
sampler = APMCABC([ff], [distance_calculator], backend, kernel, seed = 1)
steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file =4, 100, 1, 0.1, 0.03, 2, 1.0, None
print('APMCABC Inferring')
# We use resultfakeobs1 as our observed dataset
journal_apmcabc = sampler.sample([resultfakeobs1], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file)
print(journal_apmcabc.posterior_mean())
journal_apmcabc.save('apmcabc_fakeobs1.jrnl')
