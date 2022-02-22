# import logging
# logging.basicConfig(level=logging.DEBUG)
import numpy as np
from abcpy.continuousmodels import Uniform
from abcpy.output import Journal
from Model import PlateletDeposition
#from abcpy.statisticslearning import TripletDistanceLearning, SemiautomaticNN
from statistic import IdentityChosen
#from abcpy.statistics import NeuralEmbedding
#from posterior_mode import compute_posterior_mode
#import pylab as plt
###############################################
# ########### Read Observed data ################
from numpy import genfromtxt
my_data = genfromtxt('AllResults.csv', delimiter=',')

sample = False
compute_MLE = False
analysis = False
plot = False
pathology = True

if sample:
    from abcpy.backends import BackendMPI
    backend = BackendMPI()
    for ind in range(1,10):
        print(ind)

        # True parameter 89.0, 76.0, 2.49, 7e-3, 7.7, 6e-3, 8e-4
        noAP, noNAP, SR_x = int(my_data[0, 16]), int(my_data[0, 11]), float(my_data[0, 21])
        # The parameters considered random and we want to infer
        pAd = Uniform([[5], [150]], name='pAD')
        pAg = Uniform([[5], [150]], name='pAg')
        pT = Uniform([[0.1], [10.0]], name='pT')
        pF = Uniform([[0.1e-3], [9.0e-3]], name='pF')
        aT = Uniform([[0], [10]], name='aT')
        v_z_AP = Uniform([[1.0e-3], [9.0e-3]], name='v_z_AP')
        v_z_NAP = Uniform([[1.0e-4], [9.0e-4]], name='v_z_NAP')
        PD = PlateletDeposition([noAP, noNAP, SR_x, pAd, pAg, pT, pF, aT, v_z_AP, v_z_NAP], name='PD')
        # XObserved = np.array([0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.11600000e+05, 5.96585366e+03,
        #                       2.00000000e+01, 3.90572997e+03, 1.59328549e+01, 1.38943902e+05, 0.00000000e+00,
        #                       6.00000000e+01, 3.42727305e+03, 2.80052570e+01, 8.57585366e+04, 0.00000000e+00,
        #                       1.20000000e+02, 2.33523014e+03, 7.57715388e+01, 4.25231707e+04, 0.00000000e+00,
        #                       3.00000000e+02, 1.74166329e+02, 2.46413793e+03, 5.15975610e+03, 0.00000000e+00])
        # obsdata = [np.array(XObserved).reshape(1, -1)]

        # # Define kernel and join the defined kernels
        from abcpy.perturbationkernel import DefaultKernel
        kernel = DefaultKernel([pAd, pAg, pT, pF, aT, v_z_AP, v_z_NAP])

        # Define Distance functions
        from abcpy.distances import Euclidean
        from statistic import Multiply
        L = np.load('Data/L_all_3_cross.npz')['L']
        stat_mult = Multiply(L=L, degree=3, cross=True)
        dist_calc_mult = Euclidean(stat_mult)

        # SABC - Multiply##
        from abcpy.inferences import SABC
        print('Inference using Classifier Loss')
        sampler = SABC([PD], [dist_calc_mult], backend, kernel, seed=1)
        steps, epsilon, n_samples, n_samples_per_param, ar_cutoff, full_output, journal_file = 20, 10e20, 511, 1, 0.001, 1, None
        print('SABC Inferring')


        fakedata = PD.forward_simulate([noAP, noNAP, SR_x, 89.0, 76.0, 2.49, 7e-3, 7.7, 6e-3, 8e-4], k=1)
        print(fakedata)
        journal_sabc = sampler.sample(observations=[fakedata], steps=steps, epsilon=epsilon, n_samples=n_samples,
                                      n_samples_per_param=n_samples_per_param, beta=2, delta=0.2, v=0.3,
                                      ar_cutoff=0.001, resample=None, n_update=None, full_output=1,
                                      journal_file=journal_file)
        print(journal_sabc.posterior_mean())
        journal_sabc.save('Data/Journals/sabc_obs_fake_'+str(ind)+'_multiply.jrnl')


if compute_MLE:
    from posterior_mode import compute_posterior_mode
    joint_posterior_mode = np.zeros(shape=(10, 7))
    marginal_posterior_mode = np.zeros(shape=(10, 7))
    for whichobs in range(10):
        print(whichobs)
        joint_posterior_mode[whichobs, :] = compute_posterior_mode(whichobs, dim=None, type='multiply', data_type='_obs_fake_')
        for dim_ind in range(7):
            marginal_posterior_mode[whichobs,dim_ind] = compute_posterior_mode(whichobs, dim=dim_ind, type='multiply', data_type='_obs_fake_')
        print(joint_posterior_mode[whichobs, :], marginal_posterior_mode[whichobs, :])
        np.savez('posterior_modes_fake_multiply', joint_posterior_mode=joint_posterior_mode, marginal_posterior_mode=marginal_posterior_mode)

if analysis:
    true_value = [89.0, 76.0, 2.49, 7e-3, 7.7, 6e-3, 8e-4]
    joint_posterior_mode = np.load('Data/posterior_modes_fake_multiply.npz')['joint_posterior_mode']
    print(joint_posterior_mode)
    for ind in range(7):
        print(np.mean(joint_posterior_mode[:, ind])-3*np.std(joint_posterior_mode[:,ind]), np.mean(joint_posterior_mode[:, ind])+3*np.std(joint_posterior_mode[:,ind]))
        print(true_value[ind])

if pathology:
    ### Some definitions ###
    renewed = list(range(6)) + list(range(23, 65))
    mode = np.load('Data/posterior_mode_multiply.npz')['joint_posterior_mode'][renewed,:]
    print((np.std(mode[:16,:],axis=0)[1]+np.std(mode[16:32, :], axis=0)[1]+np.std(mode[32:, :], axis=0)[1])/3)
    print((np.std(mode[:16, :], axis=0)[4] + np.std(mode[16:32, :], axis=0)[4] + np.std(mode[32:, :], axis=0)[4]) / 3)
    print((np.std(mode[:16, :], axis=0)[6] + np.std(mode[16:32, :], axis=0)[6] + np.std(mode[32:, :], axis=0)[6]) / 3)
    #print(mode.shape)
    first_group_median = np.median(mode[:16,:],axis=0)
    second_group_median = np.median(mode[16:32,:],axis=0)
    third_group_median = np.median(mode[32:,:],axis=0)
    #print(first_group_median, second_group_median, third_group_median)
    distances_eu = np.zeros(shape=(48, 7, 3))
    distances_abs = np.zeros(shape=(48, 7, 3))
    for ind in range(48):
        distances_eu[ind, :, 0] = np.power(mode[ind,:] - first_group_median, 2)
        distances_eu[ind, :, 1] = np.power(mode[ind, :] - second_group_median, 2)
        distances_eu[ind, :, 2] = np.power(mode[ind, :] - third_group_median, 2)
        distances_abs[ind, :, 0] = np.abs(mode[ind,:] - first_group_median)
        distances_abs[ind, :, 1] = np.abs(mode[ind, :] - second_group_median)
        distances_abs[ind, :, 2] = np.abs(mode[ind, :] - third_group_median)
    # Between healthy vs dialysis
    FP1, TP1, FN1, TN1 = 0, 0, 0, 0
    for ind in range(16,48):
        #print(ind)
        if distances_abs[ind, 4, 1] >= distances_abs[ind, 4, 2]:
            #print('belongs to healthy')
            if ind <32:
                FN1 += 1
            else:
                TN1 += 1
        else:
            #print('belongs to dialysis')
            if ind <32:
                TP1 += 1
            else:
                FP1 += 1
    sens = TP1/(TP1+FN1)
    print('Sensitivity :'+ str(sens))
    spec = TN1 / (TN1 + FP1)
    print('Specificity :' + str(spec))
    print('Positive Predictive Value :' + str(TP1/(TP1 + FP1)))
    print('Negative Predictive Value :' + str(TN1 / (TN1 + FN1)))
    print('Positive LHD ratio :'+str(sens / (1-spec)))
    print('Negative LHD ratio :' + str((1 - sens) / (spec)))
    # Between healthy vs COPD
    FP2, TP2, FN2, TN2 = 0, 0, 0, 0
    for ind in list(range(16)) + list(range(32, 48)):
        #print(ind)
        if distances_abs[ind, 1, 0] >= distances_abs[ind, 1, 2]:
            #print('belongs to healthy')
            if ind <16:
                FN2 += 1
            else:
                TN2 += 1
        else:
            #print('belongs to copd')
            if ind <16:
                TP2 += 1
            else:
                FP2 += 1
    print(FP2, TP2, FN2, TN2)
    sens = TP2 / (TP2 + FN2)
    print('Sensitivity :' + str(sens))
    spec = TN2 / (TN2 + FP2)
    print('Specificity :' + str(spec))
    print('Positive Predictive Value :' + str(TP2 / (TP2 + FP2)))
    print('Negative Predictive Value :' + str(TN2 / (TN2 + FN2)))
    print('Positive LHD ratio :' + str(sens / (1 - spec)))
    print('Negative LHD ratio :' + str((1 - sens) / (spec)))