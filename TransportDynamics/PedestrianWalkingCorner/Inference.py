import numpy as np
from abcpy.continuousmodels import Uniform
from Model import FFSimulator1, FFSimulator2, FFSimulator3

# Define the most important options, i.e. which ABC algorithm to use and which model should be evaluated
# This is make to simplify the switching operation between different models
# exp_dataset is simply the name of the output journal file containing results
abc_method = 'apmcabc'
sim_model = 'ff3'
exp_dataset = 'fakeobs1'

# Define Graphical Model
# The distribution for all parameters used in the 3 models are provided
# Not all models use all parameters, only kS, kD, kW, decay and diffusion are common to all models
# BufferRatio and BufferAngle are specific of the model by Li et al.
# A, B, P1 and P2 are specific of the model by Dias et al.
kS = Uniform([[0], [10]], name='kS')
kD = Uniform([[0], [10]], name='kD')
kW = Uniform([[0], [10]], name='kW')
decay = Uniform([[0], [1]], name='decay')
diffusion = Uniform([[0], [1]], name='diffusion')
liBufferRatio = Uniform([[0], [3]], name='liBufferRatio')
liBufferAngle = Uniform([[45], [90]], name='liBufferAngle')
diasA = Uniform([[0],[3]], name='diasA')
diasB = Uniform([[0],[3]], name='diasB')
diasP1 = Uniform([[0],[1]], name='diasP1')
diasP2 = Uniform([[0],[1]], name='diasP2')

# Define simulator and perturbation kernels
# Since each model use a different number of parameters different ifs are employed
# The same name is used for the simulator so that there is no need to distinguish hereafter
from abcpy.perturbationkernel import MultivariateNormalKernel, JointPerturbationKernel
if sim_model=='ff1':
    ff = FFSimulator1([kS, kD, kW, decay, diffusion], name = 'ff')
    kernelcontinuous = MultivariateNormalKernel([kS, kD, kW, decay, diffusion])
if sim_model=='ff2':
    ff = FFSimulator2([kS, kD, kW, decay, diffusion, liBufferRatio, liBufferAngle], name = 'ff')
    kernelcontinuous = MultivariateNormalKernel([kS, kD, kW, decay, diffusion, liBufferRatio, liBufferAngle])
if sim_model=='ff3':
    ff = FFSimulator3([kS, kD, kW, decay, diffusion, diasA, diasB, diasP1, diasP2], name = 'ff')
    kernelcontinuous = MultivariateNormalKernel([kS, kD, kW, decay, diffusion, diasA, diasB, diasP1, diasP2])
kernel = JointPerturbationKernel([kernelcontinuous])

# Example to Generate Data to check forward simulations work without problems
if sim_model=='ff1':
    fftry = FFSimulator1([4.0, 1.0, 1.5, 0.1, 0.25], name = 'fftry')
    resultfakeobs1 = fftry.forward_simulate([4.0, 1.0, 1.5, 0.1, 0.25], 10)
    resultfakeobs2 = fftry.forward_simulate([5.0, 1.0, 1.5, 0.1, 0.25], 10)
if sim_model=='ff2':
    fftry = FFSimulator2([4.0, 1.0, 1.5, 0.1, 0.25, 0.5, 50.0], name = 'fftry')
    resultfakeobs1 = fftry.forward_simulate([4.0, 1.0, 1.5, 0.1, 0.25, 0.5, 50.0], 10)
    resultfakeobs2 = fftry.forward_simulate([5.0, 1.0, 1.5, 0.1, 0.25, 0.5, 50.0], 10)
if sim_model=='ff3':
    fftry = FFSimulator3([4.0, 1.0, 1.5, 0.1, 0.25, 1.0, 2.0, 0.725, 0.349], name = 'fftry')
    resultfakeobs1 = fftry.forward_simulate([4.0, 1.0, 1.5, 0.1, 0.25, 1.0, 2.0, 0.725, 0.349], 10)
    resultfakeobs2 = fftry.forward_simulate([5.0, 1.0, 1.5, 0.1, 0.25, 1.0, 2.0, 0.725, 0.349], 10)

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
# DistanceType1 is the Euclidean distance between the heatmaps
# DistanceType2 is the Euclidean distance between the sorted heatmaps
# DistanceType3 is a measure which determine how close are the peaks both in maximum size and in location
# DistanceType4 is the Euclidean distance between pedestrian position at each time step
from Distance import DistanceType1
distance_calculator = DistanceType1(statistics_calculator)

# Check whether the distance works
if distance_calculator.distance(resultfakeobs1, resultfakeobs1)==distance_calculator.distance(resultfakeobs1, resultfakeobs2):
    print('Something may be wrong with the distance!')

###############################################################################
#                                APMCABC                                      #
###############################################################################

if abc_method=='apmcabc':
    from abcpy.inferences import APMCABC
    sampler = APMCABC([ff], [distance_calculator], backend, kernel, seed = 1)
    steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file = 4, 1000, 1, 0.1, 0.03, 2, 1.0, None
    print('APMCABC Inferring')
    
    # We use resultfakeobs1 as our observed dataset
    journal_apmcabc = sampler.sample([resultfakeobs1], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file)
    print(journal_apmcabc.posterior_mean())
    journal_apmcabc.save('apmcabc_' + sim_model + '_' + exp_dataset + '.jrnl')

###############################################################################
#                                  SABC                                       #
###############################################################################
if abc_method=='sabc':
    from abcpy.inferences import SABC
    sampler = SABC([ff], [distance_calculator], backend, kernel, seed = 1)
    steps, epsilon, n_samples, n_samples_per_param, beta, delta, v, ar_cutoff, resample, n_update, adaptcov, full_output, journal_file = 2, 40, 1000, 1, 2, 0.2, 0.3, 0.5, None, None, 1, 0, None
    print('SABC Inferring')
    
    ## We use resultfakeobs1 as our observed dataset
    journal_sabc1 = sampler.sample([resultfakeobs1], steps, epsilon, n_samples, n_samples_per_param, beta, delta, v, ar_cutoff, resample, n_update, adaptcov, full_output, journal_file)
    print(journal_sabc1.posterior_mean())
    journal_sabc1.save('sabc_' + sim_model + '_' + exp_dataset + '.jrnl')
