import numpy as np
import time
from Model import Helium

# Define Grphical Model
from abcpy.continuousmodels import Uniform
sigma = Uniform([[0.1], [0.8]], name='sigma')
epsilon = Uniform([[0.01], [1.0]], name='epsilon')
relentropy = Helium([sigma, epsilon], name = 'relentropy')

# Simulate Dataset
result = relentropy.forward_simulate([0.2556, 0.141], 1)[0]
np.save('Result/obs_data_'+str(ind), result)

# Define backend
from abcpy.backends import BackendMPI as Backend
backend = Backend()

# Define Statistics
from abcpy.statistics import Identity
statistics_calculator = Identity(degree=1, cross=False)

# Define distance
from KLdistance import KLdistance
distance_calculator = KLdistance(statistics_calculator)

# Define kernel
from abcpy.perturbationkernel import DefaultKernel
kernel = DefaultKernel([sigma, epsilon])


## APMCABC ##
from abcpy.inferences import APMCABC
sampler = APMCABC([relentropy], [distance_calculator], backend, kernel, seed = 1)
steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file =10, 1000, 1, 0.1, 0.03, 2.0, 1, None

# Import Simulated Dataset
relentropy_obs = [np.load('Result/obs_data.npy')]
print('APMCABC Inferring Heliuma Potential')
journal_apmcabc = sampler.sample([relentropy_obs], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file)
journal_apmcabc.save('Result/APMCABCHelium.jrnl')

