import numpy as np
from abcpy.statistics import Identity
from Distance import Abserror
from abcpy.inferences import APMCABC
from abcpy.continuousmodels import Uniform
from Model import TIP4PGromacsOOOH as Water

# Define Graphical model
theta1 = Uniform([[.281] , [.53]],name='theta1')
theta2 = Uniform([[0.2] , [0.9]],name='theta2')
water = Water([theta1, theta2, 2500000])

# Define distance and statistics
statistics_calculator = Identity(degree = 1, cross = False)
distance_calculator = Abserror(statistics_calculator)

# Define kernel
from abcpy.backends import BackendMPI as Backend
backend = Backend()
from abcpy.perturbationkernel import DefaultKernel
kernel = DefaultKernel([theta1, theta2])


######### Inference for simulated data ###############
water_obs = [np.load('Data/obs_data.npy')]

sampler = APMCABC([water], distance_calculator, backend, kernel, seed = 1)
steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file =10, 100, 1, 0.1, 0.03, 2.0, 1, None

print('TIP4P: APMCABC Inferring for simulated data')
journal_apmcabc = sampler.sample([water_obs], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file)
print('TIP4P: APMCABC done for simulated data')
journal_apmcabc.save('Result/MD_GROMACS_APMCABC_obs.jrnl')

######### Inference for Experimental data 1 (Neutron Diffraction of Water) ###############
water_obs = [np.load('Data/exp_data.npy')]

sampler = APMCABC([water], distance_calculator, backend, kernel, seed = 1)
steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file =10, 100, 1, 0.1, 0.03, 2.0, 1, None

print('TIP4P: APMCABC Inferring for experimental data 1')
journal_apmcabc = sampler.sample([water_obs], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file)
print('TIP4P: APMCABC done for experimental data 1')
journal_apmcabc.save('Result/MD_GROMACS_APMCABC_exp.jrnl')


######### Inference for Experimental data 2 (X-ray Diffraction of Water) ###############
from Model import TIP4PGromacsOO as WaterOO
theta1 = Uniform([[.281] , [.53]],name='theta1')
theta2 = Uniform([[0.2] , [0.9]],name='theta2')

# Define Graphical model
water = WaterOO([theta1, theta2, 2500000])

water_obs = [np.load('Data/exp_data_2.npy')]

sampler = APMCABC([water], distance_calculator, backend, kernel, seed = 1)
steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file =10, 100, 1, 0.1, 0.03, 2.0, 1, None

print('TIP4P: APMCABC Inferring for experimental data 2')
journal_apmcabc = sampler.sample([water_obs], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file)
print('TIP4P: APMCABC done for experimental data 2')
journal_apmcabc.save('Result/MD_GROMACS_APMCABC_exp_2.jrnl')

