# import logging
# logging.basicConfig(level=logging.DEBUG)
import numpy as np
from abcpy.continuousmodels import Uniform
from abcpy.output import Journal
from Model import PlateletDeposition
from abcpy.statisticslearning import TripletDistanceLearning, SemiautomaticNN
from statistic import IdentityChosen
from abcpy.statistics import NeuralEmbedding
from posterior_mode import compute_posterior_mode
import pylab as plt
###############################################
# ########### Read Observed data ################
from numpy import genfromtxt
my_data = genfromtxt('Data/AllResults.csv', delimiter=',')
fake = False
train_classifier = False
simulate_pilot = False
train_save_NN = False
sample = False
plot = True
predict = False
###############################################
# Define backend
from abcpy.backends import BackendMPI, BackendDummy
backend = BackendDummy()
#backend = BackendMPI()
print('Hello')
renewed = list(range(6)) + list(range(23,65))
if train_classifier:
    alldata, label, infoindi = [], [], []
    for whichobs in renewed:
        XObserved = [np.hstack((np.array([0, 20, 60, 120, 300]).reshape(5, 1),
                                my_data[:, 1:21][whichobs, :].reshape(4, 5).transpose())).flatten().reshape(1, -1)][0]
        alldata.append(XObserved[0])
        obsdata = [np.array(XObserved).reshape(1, -1)]
        # Define the indices used as features to remove timestamp etc. and values with NA
        AllIndices = list(np.arange(0, 25, 1))
        NANTimeIndices = list(np.argwhere(np.isnan(XObserved[0, :])).reshape(-1, )) + [0, 5, 10, 15, 20] + list(np.arange(0, 5, 1)) + [4, 9, 14, 19, 24]
        InformativeIndices = [item for item in AllIndices if item not in NANTimeIndices]
        infoindi.append(InformativeIndices)
        noAP, noNAP, SR_x = int(my_data[whichobs, 16]), int(my_data[whichobs, 11]), float(my_data[whichobs, 21])
        label.append(my_data[whichobs, 23])
    result = set(infoindi[0]).intersection(*infoindi[1:])
    print(result)
    alldata_cleaned = np.array(alldata)[:,list(result)]
    from abcpy.statistics import Identity
    stat = Identity(degree=3, cross=True)
    XChosen = stat.statistics([[alldata_cleaned[i,:]] for i in range(alldata_cleaned.shape[0])])
    print(XChosen.shape)
    from metric_learn import LMNN
    metric = LMNN(init='auto',k=6, min_iter=10000, max_iter=50000, convergence_tol=1e-6, learn_rate=1e-10, regularization=.5, n_components = 2)
    metric.fit(XChosen, label)
    L = metric.components_
    np.savez('Data/L_all_3_cross.npz', L=L)

    L = np.load('Data/L_all_3_cross.npz')['L']
    X_lmnn = XChosen.dot(L.T)
    print(X_lmnn.shape)
    import pylab as plt
    plt.figure()
    plt.plot(X_lmnn[:16,0], X_lmnn[:16,1], 'k*', label='Patients with COPD')
    plt.plot(X_lmnn[16:32,0], X_lmnn[16:32,1], 'r*', label='Patients having dialysis')
    plt.plot(X_lmnn[32:,0], X_lmnn[32:,1], 'b*', label ='Healthy volunteers')
    plt.xlabel(r'$X_1$',fontsize=20)
    plt.ylabel(r'$X_2$',fontsize=20)
    #plt.title('2-dim projected space', fontsize=20)
    plt.legend(fontsize=10)
    plt.savefig('FinalFigures/Fig6.pdf')

    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics.cluster import adjusted_rand_score
    cluster_lmnn = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
    cluster_lmnn.fit_predict(X_lmnn)
    print('AR score of Euclidean distance btwn LMNN: '+str(adjusted_rand_score(label, cluster_lmnn.labels_)))

for whichobs in renewed:
    ########### Define Graphical Model ############
    #Define which experimental study
    if fake:
        # True parameter 89.0, 76.0, 2.49, 7e-3, 7.7, 6e-3, 8e-4
        XObserved = np.array([0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.11600000e+05, 5.96585366e+03,
         2.00000000e+01, 3.90572997e+03, 1.59328549e+01, 1.38943902e+05, 0.00000000e+00,
         6.00000000e+01, 3.42727305e+03, 2.80052570e+01, 8.57585366e+04, 0.00000000e+00,
         1.20000000e+02, 2.33523014e+03, 7.57715388e+01, 4.25231707e+04, 0.00000000e+00,
         3.00000000e+02, 1.74166329e+02, 2.46413793e+03, 5.15975610e+03, 0.00000000e+00])
        obsdata = [np.array(XObserved).reshape(1, -1)]
        noAP, noNAP, SR_x = int(my_data[0, 16]), int(my_data[0, 11]), float(my_data[0, 21])
        # Define the indices used as features to remove timestamp etc. and values with NA
        AllIndices = list(np.arange(0,25,1))
        NANTimeIndices = [0, 5, 10, 15, 20] + list(np.arange(0, 5, 1)) + [4, 9, 14, 19, 24]
        InformativeIndices = [item for item in AllIndices if item not in NANTimeIndices]
    else:
        XObserved = [np.hstack((np.array([0,20,60,120,300]).reshape(5,1),my_data[:,1:21][whichobs,:].reshape(4, 5).transpose())).flatten().reshape(1, -1)][0]
        obsdata = [np.array(XObserved).reshape(1, -1)]
        # Define the indices used as features to remove timestamp etc. and values with NA
        AllIndices = list(np.arange(0, 25, 1))
        NANTimeIndices = list(np.argwhere(np.isnan(XObserved[0, :])).reshape(-1, )) + [0, 5, 10, 15, 20]
        InformativeIndices = [item for item in AllIndices if item not in NANTimeIndices]
        noAP, noNAP, SR_x = int(my_data[whichobs, 16]), int(my_data[whichobs, 11]), float(my_data[whichobs, 21])

    ## Summary statistics taking informative indices
    identity = IdentityChosen(InformativeIndices=InformativeIndices, degree=1, cross=False)  # to apply before computing the statistics
    ## Random values
    # The parameters considered random and we want to infer
    pAd = Uniform([[5], [150]], name='pAD')
    pAg = Uniform([[5], [150]], name='pAg')
    pT = Uniform([[0.1], [10.0]], name='pT')
    pF = Uniform([[0.1e-3], [9.0e-3]], name='pF')
    aT = Uniform([[0], [10]], name='aT')
    v_z_AP =  Uniform([[1.0e-3], [9.0e-3]], name='v_z_AP')
    v_z_NAP =  Uniform([[1.0e-4], [9.0e-4]], name='v_z_NAP')
    PD = PlateletDeposition([noAP, noNAP, SR_x, pAd, pAg, pT, pF, aT, v_z_AP, v_z_NAP], name = 'PD')
    #fakedata = PD.forward_simulate([noAP, noNAP, SR_x, 89.0, 76.0, 2.49, 7e-3, 7.7, 6e-3, 8e-4], k=1)
    # # Define kernel and join the defined kernels
    from abcpy.perturbationkernel import DefaultKernel
    kernel = DefaultKernel([pAd, pAg, pT, pF, aT, v_z_AP, v_z_NAP])

    # Now learn the optimal summary statistics to be used
    ########## Define distance based on embedding learned #####
    if simulate_pilot:
        print('Simulate')
        from Model import DrawFromPrior
        draw_from_prior = DrawFromPrior([PD], backend=backend)
        parameters, simulations = draw_from_prior.sample(255, n_samples_per_param=1)
        if fake:
            np.savez('Data/Pilots/simulation_pilot_fake.npz', parameters=parameters, simulations=simulations)
        else:
            np.savez('Data/Pilots/simulation_pilot_'+str(whichobs)+'.npz', parameters=parameters, simulations=simulations)

    if train_save_NN:
        # # Define backend
        # from abcpy.backends import BackendDummy
        # backend = BackendDummy()
        if fake:
            parameters = np.load('Data/Pilots/simulation_pilot_fake.npz')['parameters']
            simulations = np.load('Data/Pilots/simulation_pilot_fake.npz')['simulations']
        else:
            parameters = np.load('Data/ailots/simulation_pilot_'+str(whichobs)+'.npz')['parameters']
            simulations = np.load('Data/Pilots/simulation_pilot_'+str(whichobs)+'.npz')['simulations']
        print(parameters.shape, simulations.shape)
        #now train the NNs with the different methods with the generated data
        print("semiNN")
        semiNN = SemiautomaticNN([PD], identity, backend=backend, parameters=parameters, simulations=simulations,
                                  early_stopping=False,  batch_size=32,# early stopping
                                  seed=1, n_epochs=1000, lr=1e-4, scale_samples=False) #1000
        print("triplet")
        triplet = TripletDistanceLearning([PD], identity, backend=backend, parameters=parameters, simulations=simulations,
                                           early_stopping=False, batch_size=32, # early stopping
                                           seed=1, n_epochs=2000, lr=1e-3, scale_samples=False) # 2000
        # get the statistics from the StatisticsLearning object:
        learned_seminn_stat = semiNN.get_statistics()
        learned_triplet_stat = triplet.get_statistics()
        # saving the learned net
        if fake:
            learned_seminn_stat.save_net("Data/Pilots/seminn_net_fake.pth")
            learned_triplet_stat.save_net("Data/Pilots/triplet_net_fake.pth")
        else:
            learned_seminn_stat.save_net("Data/Pilots/seminn_net_"+str(whichobs)+".pth")
            learned_triplet_stat.save_net("Data/Pilots/triplet_net_"+str(whichobs)+".pth")

    if sample:
        if fake:
            learned_seminn_stat_loaded = NeuralEmbedding.fromFile("Data/Pilots/seminn_net_fake.pth", input_size=len(InformativeIndices), output_size=7, previous_statistics=identity)
            learned_triplet_stat_loaded = NeuralEmbedding.fromFile("Data/Pilots/triplet_net_fake.pth", input_size=len(InformativeIndices), output_size=7, previous_statistics=identity)
        else:
            learned_seminn_stat_loaded = NeuralEmbedding.fromFile("Data/Pilots/seminn_net_"+str(whichobs)+".pth", input_size=len(InformativeIndices), output_size=7, previous_statistics=identity)
            learned_triplet_stat_loaded = NeuralEmbedding.fromFile("Data/Pilots/triplet_net_"+str(whichobs)+".pth", input_size=len(InformativeIndices), output_size=7, previous_statistics=identity)

        # Define Distance functions
        from abcpy.distances import Euclidean
        dist_calc_seminn = Euclidean(learned_seminn_stat_loaded)
        dist_calc_triplet = Euclidean(learned_triplet_stat_loaded)
        from statistic import Multiply
        L = np.load('Data/L_all_3_cross.npz')['L']
        stat_mult = Multiply(L=L, degree=3, cross=True)
        dist_calc_mult = Euclidean(stat_mult)
        print(dist_calc_mult.distance(obsdata, obsdata))

        # print('Inference starting')
        # # ## SABC - SemiNN##
        from abcpy.inferences import SABC
        print('Inference using Semi NN')
        sampler = SABC([PD], [dist_calc_seminn], backend, kernel, seed=1)
        steps, epsilon, n_samples, n_samples_per_param, ar_cutoff, full_output, journal_file = 25, 10e20, 511, 1, 0.001, 1, None
        print('SABC Inferring')
        journal_sabc = sampler.sample(observations=[obsdata], steps=steps, epsilon=epsilon, n_samples=n_samples,
                                       n_samples_per_param=n_samples_per_param, beta=2, delta=0.2, v=0.3,
                                       ar_cutoff=0.001, resample=None, n_update=None, full_output=1,
                                       journal_file=journal_file)
        print(journal_sabc.posterior_mean())
        if fake:
            journal_sabc.save('Data/Journals/sabc_obs_fake_seminn.jrnl')
        else:
            journal_sabc.save('Data/Journals/sabc_obs_' + str(whichobs) + '_seminn.jrnl')
        ## SABC - Triplet##
        print('Inference using Triplet Loss')
        sampler = SABC([PD], [dist_calc_triplet], backend, kernel, seed=1)
        steps, epsilon, n_samples, n_samples_per_param, ar_cutoff, full_output, journal_file = 10, 10e20, 511, 1, 0.001, 1, None
        print('SABC Inferring')
        journal_sabc = sampler.sample(observations=[obsdata], steps=steps, epsilon=epsilon, n_samples=n_samples,
                                       n_samples_per_param=n_samples_per_param, beta=2, delta=0.2, v=0.3,
                                       ar_cutoff=0.001, resample=None, n_update=None, full_output=1,
                                       journal_file=journal_file)
        print(journal_sabc.posterior_mean())
        if fake:
            journal_sabc.save('Data/Journals/sabc_obs_fake_triplet.jrnl')

        else:
            journal_sabc.save('Data/Journals/sabc_obs_' + str(whichobs) + '_triplet.jrnl')
        # SABC - Multiply##
        print('Inference using Classifier Loss')
        sampler = SABC([PD], [dist_calc_mult], backend, kernel, seed=1)
        steps, epsilon, n_samples, n_samples_per_param, ar_cutoff, full_output, journal_file = 20, 10e20, 511, 1, 0.001, 1, None
        print('SABC Inferring')
        journal_sabc = sampler.sample(observations=[obsdata], steps=steps, epsilon=epsilon, n_samples=n_samples,
                                      n_samples_per_param=n_samples_per_param, beta=2, delta=0.2, v=0.3,
                                      ar_cutoff=0.001, resample=None, n_update=None, full_output=1,
                                      journal_file=journal_file)
        print(journal_sabc.posterior_mean())
        if fake:
            journal_sabc.save('Data/Journals/sabc_obs_fake_multiply.jrnl')
        else:
            journal_sabc.save('Data/Journals/sabc_obs_' + str(whichobs) + '_multiply.jrnl')

    if plot:
        x_val = [20, 120, 300]
        InformativeIndices = [[6, 16, 21], [7, 17, 22], [8, 18, 23]]
        nameslabel = [r'$\mathbb{N}_{agg-clust} $', r'$\mathcal{S}_{agg-clust}}$', r'$\mathbb{N}_{platelet}$']
        filename = ['naggclust', 'saggclust', 'np']
        if fake:
            for ind in range(3):
                fig, ax = plt.subplots()
                ax.set_title(nameslabel[ind], fontsize=35)
                ax.plot(x_val, obsdata[0][0, InformativeIndices[ind]], color="black", marker=".",
                        label=r'$\mathcal{x}^0$')
                ax.set_xlabel("Time " + r"$(s)$", fontsize=30)
                ax.set_xticks(x_val)
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                plt.xticks(fontsize=30)
                plt.yticks(fontsize=30)
            plt.tight_layout()
            plt.savefig('FinalFigures/original_data/original_data_fake'+filename[ind]+'.pdf')
            plt.close()
        else:
            fig, axs = plt.subplots(1, 3, figsize=(30, 10))  # 1
            for ind, ax in enumerate(axs.ravel()):  # 2
                ax.set_title(nameslabel[ind], fontsize=35)
                ax.plot(x_val, obsdata[0][0, InformativeIndices[ind]], color="black", marker=".",
                        label=r'$\mathcal{x}^0$')
                ax.set_xlabel("Time " + r"$(s)$", fontsize=30)
                ax.set_xticks(x_val)
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                ax.tick_params(axis='both', labelsize=30)
            plt.tight_layout()
            plt.savefig('FinalFigures/original_data/original_data_' + str(whichobs) + '.pdf')
            plt.close()