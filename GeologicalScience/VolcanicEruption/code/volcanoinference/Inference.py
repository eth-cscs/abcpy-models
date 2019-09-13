import numpy as np
from abcpy.continuousmodels import Uniform
from Model import AshDispersal

u0 = Uniform([[100],[300]], name='u0')
l0 = Uniform([[30],[100]], name='l0')
AD = AshDispersal([u0, l0], name='AD')

from abcpy.backends import BackendDummy as Backend
backend = Backend(process_per_model=2)

from abcpy.perturbationkernel import DefaultKernel
kernel = DefaultKernel([u0, l0])

from Distance import DepositionDistance
distance_calculator = DepositionDistance()

print(distance_calculator.distance("test.h5", "test.h5"))

from abcpy.inferences import SABC
sampler = SABC([AD], [distance_calculator], backend, kernel, seed = 1)
#steps, epsilon, n_samples, n_samples_per_param, beta, delta, v, ar_cutoff, resample, n_update, adaptcov, full_output = 3, [50], 50, 1, 2, 0.2, 0.3, 0.1, None, None, 1, 1
steps, epsilon, n_samples, n_samples_per_param, beta, delta, v, ar_cutoff, resample, n_update, adaptcov, full_output = 3, [1], 1, 1, 2, 0.2, 0.3, 0.1, None, None, 1, 1
print('SABC Inferring')
journal_sabc = sampler.sample("test.h5", steps, epsilon, n_samples, n_samples_per_param, beta, delta, v, ar_cutoff, resample, n_update, adaptcov, full_output)
print('SABC done')