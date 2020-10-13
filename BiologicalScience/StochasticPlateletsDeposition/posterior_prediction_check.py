import logging, os
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

L = np.array([[ 2.61836732e-05, 1.24366184e-04, 1.65418571e-04, -2.76891651e-01, -3.16512105e-05, 1.31921975e-04, 1.17853188e-02, 1.07568061e-03, -1.48116765e-02, -1.39082174e-05],
 [-5.95664564e-06, 1.97248981e-05, -1.17200781e-04, -2.66547093e-01, 6.23778123e-06, -1.93221396e-03, -4.88895269e-02, 1.80570575e-03, 8.92376774e-02, -3.33985710e-06]])
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


filename = 'fakeobs_' + str(whichobsinfer) +'.npz'
if os.path.isfile(filename):
    print('Hello')
    resultfakeobs = np.load(filename)['arr_0']
    print(len(resultfakeobs))
else:
    filenamere = '../APMCABC_Results/apmcabc_obs_' + str(whichobsinfer) + '_reweighted.jrnl'
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
        resultfakeobs.append(
            PD.forward_simulate([noAP, noNAP, SR_x, data[index_to_simulate[0], 0], data[index_to_simulate[0], 1],
                                 data[index_to_simulate[0], 2], data[index_to_simulate[0], 3],
                                 data[index_to_simulate[0], 4],
                                 data[index_to_simulate[0], 5], data[index_to_simulate[0], 6]], 1))
    print(resultfakeobs)
    np.savez('fakeobs_34', resultfakeobs)


posterior_prediction = np.zeros(shape=(100,25))
for ind in range(100):
    posterior_prediction[ind, :] = resultfakeobs[ind].squeeze().reshape(1,25)

observed = obsdata[0][:,indices_chosen]
posterior_prediction = posterior_prediction[:, indices_chosen]

print(observed.shape)
print(posterior_prediction.shape)

import pylab as plt
plt.figure()
plt.plot(np.arange(1, 11), observed.T, 'k.-')
plt.boxplot(posterior_prediction)
# plt.plot(np.max(posterior_prediction.T, axis=1))
# plt.plot(np.median(posterior_prediction.T, axis=1))
# plt.plot(np.min(posterior_prediction.T, axis=1))
plt.show()