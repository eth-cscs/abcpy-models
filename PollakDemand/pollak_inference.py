import numpy as np
from abcpy.backends import BackendDummy as Backend
from abcpy.continuousmodels import Uniform
from abcpy.inferences import APMCABC, SABC
from abcpy.perturbationkernel import (DefaultKernel, JointPerturbationKernel,
                                      MultivariateStudentTKernel)
from abcpy.statistics import Identity
from numpy.random import RandomState

import distances
from model import PollakDemand
random_state =25 

# Define Graphical Model
# TODO: get better ranges for the lower and upper bounds on these priors. Currently just random values

gamma = Uniform([[0.1],[10]], name = 'gamma') #gamma always less than output x
delta = Uniform([[0.1],[10]], name = 'delta')
sigma = Uniform([[0.1],[10]], name = 'sigma' )
ln_loc = Uniform( [[10],[15]], name = 'ln_loc')
ln_scale = Uniform([[2],[8]], name = 'ln_scale')
ln_scale_underlying = Uniform([[1],[10]], name='ln_scale_underlying')


pollak = PollakDemand([gamma, delta, sigma, ln_loc, ln_scale, ln_scale_underlying], name='Pollak', markup_or_return='return')
#pollak = PollakDemand([gamma, delta, sigma, ln_loc, ln_scale, ln_scale_underlying], name='Pollak', markup_or_return='markup')

# Define backend
backend = Backend()

#Define statistics calculator
statistics_calculator = Identity(degree = 0, cross=False) #Defined but not used in distances that do not use summary statistic

#Define Distance - Choose 1 of 3 below
distance_calculator = distances.KullbackLiebler(statistics_calculator)
#distance_calculator = distances.Hellinger(statistics_calculator)
#distance_calculator = distances.Wasserstein_distance(statistics_calculator)

# Define Kernel
kernel = DefaultKernel( [gamma, delta, sigma,ln_scale, ln_loc, ln_scale_underlying] )
# kernel_pollak_params = MultivariateStudentTKernel([gamma, delta, sigma]) 
# kernel_lognormal_params = MultivariateStudentTKernel([ln_loc, ln_scale, ln_scale_underlying])
# kernel = JointPerturbationKernel([kernel_pollak_params, kernel_lognormal_params])

# Define fake observations
fake_output_data = pollak.forward_simulate([0.5, 3, 3.1, 8, 12.0, 2.1])

## SABC ##
sampler_sabc =  SABC([pollak], [distance_calculator], backend, kernel, seed=random_state )
steps, epsilon, n_samples, n_samples_per_param, beta, delta, v, ar_cutoff, resample, n_update, full_output =2, 10.0, 4, 1, 2, 0.2, 0.3, 0.1, None, None, 1
# #adaptcov = 1
print('SABC Inferring')
journal_sabc = sampler_sabc.sample( [fake_output_data], steps, epsilon, n_samples, n_samples_per_param, beta, delta, v, ar_cutoff, resample, n_update, full_output)
print('SABC done')
journal_sabc.save('sabc_fake_data.jrnl')

## APMCABC ##
# sampler_apmcabc = APMCABC([pollak], [distance_calculator], backend, kernel, seed = random_state)
# steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file =2, 8, 1, 0.25, 0.03, 2, 1.0, None
# print('APMCABC Inferring')
# journal_apmcabc = sampler_apmcabc.sample([fake_output_data], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file) 
# journal_apmcabc.save('apmcabc_fake_data.jrnl')
