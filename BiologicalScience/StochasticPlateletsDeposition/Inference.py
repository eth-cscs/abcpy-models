import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
from abcpy.continuousmodels import Uniform
from abcpy.output import Journal
from Model import PlateletDeposition
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score
# ###############################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("whichobsinfer", help="which observation", type=int)
args = parser.parse_args()
whichobsinfer = args.whichobsinfer
# ########### Read Observed data ################
from numpy import genfromtxt
from metric_learn import LMNN
my_data = genfromtxt('DataCleanedObserved.csv', delimiter=',')
# The missing column has been imputed and actiavted platelets values computed using the formula
ShR = my_data[2:, 0]
X = my_data[2:,2:]

# ############ Learn distance from the dataset#####
# # Consider 3 classes as volunteer (16 of them), Patient Dialysé (non Diabètiqu + Diabètiqu) (16 of them) and  Patient BPCO (stable + Exacerbé) (8 of them)
# n_clus = 3
A, B, C = X[:16, :], X[16:32, :], X[32:, :]
X = np.vstack((A, B, C))
AI, BI, CI = np.ones(shape=(16,1)), np.ones(shape=(16,1)), np.ones(shape=(8,1))
Y = np.vstack((AI, 2*BI, 3*CI)).reshape(-1,)

XObserved = np.zeros(shape=(40,25))
for whichobs in range(40):
    XObserved[whichobs,:] = [np.hstack((np.array([0,20,60,120,300]).reshape(5,1),X[whichobs,:].reshape(4,5).transpose())).flatten().reshape(1,-1)][0]
# Define the indices used as features to remove timestamp etc. and values with NA
indices_chosen = [3, 4, 6, 7, 8, 16, 17, 21, 22, 23]
XChosen = XObserved[:,indices_chosen]

# Learn the embedding
# lmnn = LMNN(init='lda',k=6, min_iter=10000, max_iter=50000, convergence_tol=1e-6, learn_rate=1e-10, regularization=.5, n_components = 2)
# lmnn.fit(XChosen, Y)
# L = lmnn.components_

L = np.array([[ 2.61836732e-05, 1.24366184e-04, 1.65418571e-04, -2.76891651e-01, -3.16512105e-05, 1.31921975e-04, 1.17853188e-02, 1.07568061e-03, -1.48116765e-02, -1.39082174e-05],
 [-5.95664564e-06, 1.97248981e-05, -1.17200781e-04, -2.66547093e-01, 6.23778123e-06, -1.93221396e-03, -4.88895269e-02, 1.80570575e-03, 8.92376774e-02, -3.33985710e-06]])
# Report classification capabiliy of embedding
#X_lmnn = lmnn.transform(XChosen)

# import pylab as plt
# plt.figure()
# plt.plot(X_lmnn[:16,0], X_lmnn[:16,1], 'k*', label='Healthy Volunteers')
# plt.plot(X_lmnn[16:32,0], X_lmnn[16:32,1], 'r*', label='Patients having dialysis')
# plt.plot(X_lmnn[32:,0], X_lmnn[32:,1], 'b*', label ='COPD patients')
# plt.xlabel(r'$X_1$',fontsize=20)
# plt.ylabel(r'$X_2$',fontsize=20)
# plt.title('2-dim projected space', fontsize=20)
# plt.legend(fontsize=10)
# plt.savefig('Projected.eps')

# # cluster_lmnn = AgglomerativeClustering(n_clusters=n_clus, affinity='euclidean', linkage='ward')
# # cluster_lmnn.fit_predict(X_lmnn)
# # print('AR score of Euclidean distance btwn LMNN: '+str(adjusted_rand_score(Y, cluster_lmnn.labels_)))
#
# ########## Define distance based on embedding learned #####
# # Start from identity statistics
# from abcpy.statistics import Identity
# statistics_calculator = Identity(degree = 1, cross = False)
# # Redefine the statistics function
# statistics_calculator.statistics = lambda x, f1=statistics_calculator.statistics: (f1(x)[:,indices_chosen]).dot(L.T)
#
# # Define Euclidean distance on the embedding
# from abcpy.distances import Euclidean
# dist_calc = Euclidean(statistics_calculator)
# # ###############################################
# Define backend
from abcpy.backends import BackendDummy as Backend
backend = Backend()
# ########### Define Graphical Model ############
# Define which experimental study
#whichobs = 0
print(whichobsinfer)
obsdata = [np.array(XObserved[whichobsinfer,:]).reshape(1,-1)]
## Fixed values
# Initial values of AP and NAP
noAP = XObserved[whichobsinfer,4].astype(int)
noNAP = XObserved[whichobsinfer,3].astype(int)
# Fixed value of sheer rate used fror the experiment
SR_x = float(ShR[whichobsinfer]) #100
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

filenamere = '../APMCABC_Results/apmcabc_obs_'+str(whichobsinfer)+'_reweighted.jrnl'
journal = Journal.fromFile(filenamere)
parameters_to_show = ['pAD', 'pAg', 'pT', 'pF', 'aT', 'v_z_AP', 'v_z_NAP']
post_samples_dict = {name: np.concatenate(journal.get_parameters(-1)[name]) for name in parameters_to_show}
data = post_samples_dict[parameters_to_show[0]]
for ind in range(1, len(parameters_to_show)):
    data = np.hstack((data, post_samples_dict[parameters_to_show[ind]]))
weights = np.concatenate(journal.get_weights(-1))
index_to_simulate = np.random.choice(data.shape[0], 100, p=weights)

resultfakeobs = []
for ind in range(100):
    resultfakeobs.append(PD.forward_simulate([noAP, noNAP, SR_x, data[index_to_simulate[0],0], data[index_to_simulate[0],1],
                                      data[index_to_simulate[0],2], data[index_to_simulate[0],3], data[index_to_simulate[0],4],
                                      data[index_to_simulate[0],5], data[index_to_simulate[0],6]], 1))
print(resultfakeobs)
np.savez('fakeobs_34', resultfakeobs)

# #resultfakeobs1 = PD.forward_simulate([4208, 147160, 400.0, 89.23867991870651, 76.86412615447622, 2.4914170141882135, 0.0007587473496343516, 7.702218665710551, 0.006532910756608233, 0.0008366912599321585], 1)
#print(dist_calc.distance(obsdata, resultfakeobs1))


# # # Define kernel and join the defined kernels
# from abcpy.perturbationkernel import DefaultKernel
# kernel = DefaultKernel([pAd, pAg, pT, pF, aT, v_z_AP, v_z_NAP])
# #
# # ## APMCABC ##
# from abcpy.inferences import APMCABC
# sampler = APMCABC([PD], [dist_calc], backend, kernel, seed = 1)
# steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file =10, 3500, 1, 0.1, 0.001, 2, 1, None
# print('APMCABC Inferring')
# # We use resultfakeobs1 as our observed dataset
# journal_apmcabc = sampler.sample([obsdata], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file)
# print(journal_apmcabc.posterior_mean())
# journal_apmcabc.save('apmcabc_obs_'+str(whichobsinfer)+'.jrnl')
