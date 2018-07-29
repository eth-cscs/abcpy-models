import numpy as np
from abcpy.continuousmodels import Uniform
from abcpy.discretemodels import Binomial
from Model import FFSimulator, DiscreteUniform


# Define Graphical Model
bufferRatio = Uniform([[0], [10]], name='bufferratio')
bufferAngle = DiscreteUniform([[45],[90]], name='bufferAngle')
kW = Uniform([[0], [10]], name='kW')
kS = Uniform([[0], [10]], name='kS')
kD = Uniform([[0], [10]], name='kD')
decay = Uniform([[0], [1]], name='decay')
diffusion = Uniform([[0], [1]], name='diffusion')
ff = FFSimulator([bufferRatio, bufferAngle, kW, kS, kD, decay, diffusion], name = 'ff')

# Example to Generate Data to check it's correct
fftry = FFSimulator([8.2, 50, 3.2, 4.1, 1.1, .3, .1], name = 'fftry')
resultfakeobs1 = fftry.forward_simulate([8.2, 50, 3.2, 4.1, 1.1, .3, .1], 1)
resultfakeobs2 = fftry.forward_simulate([9.2, 50, 3.2, 4.1, 1.1, .3, .1], 1)

print('# Check the datasets are different')
print(np.mean(resultfakeobs1[0]))
print(np.mean(resultfakeobs2[0]))
  
# Define backend
from abcpy.backends import BackendDummy as Backend
backend = Backend()

# Define Statistics
from abcpy.statistics import Identity
statistics_calculator = Identity(degree=1, cross=False)

# Define distance
from Distance import DistanceType1
distance_calculator = DistanceType1(statistics_calculator)

print('# Check whether the distance works')
print(distance_calculator.distance(resultfakeobs1, resultfakeobs1))
print(distance_calculator.distance(resultfakeobs1, resultfakeobs2))


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
