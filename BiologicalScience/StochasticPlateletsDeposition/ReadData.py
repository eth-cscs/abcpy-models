import numpy as np
from numpy import genfromtxt
import pylab as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score
colorarray = ['tomato', 'royalblue', 'limegreen']

my_data = genfromtxt('DataCleanedObserved.csv', delimiter=',')

# The missing column has been imputed and actiavted platelets values computed using the formula
X = my_data[2:,2:]

# ############ Learn distance from the dataset#####
# # Consider 3 classes as volunteer (16 of them), Patient Dialysé (non Diabètiqu + Diabètiqu) (16 of them) and  Patient BPCO (stable + Exacerbé) (8 of them)
# n_clus = 3
A, B, C = X[:16, :], X[16:32, :], X[32:, :]
# AI, BI, CI = np.ones(shape=(16,1)), np.ones(shape=(16,1)), np.ones(shape=(8,1))
# Y = np.vstack((AI, 2*BI, 3*CI)).reshape(-1,)
X = np.vstack((A, B, C))

# Choice
indices_chosen = [3, 4, 6, 7, 8, 16, 17, 21, 22, 23]
# Start from identity statistics
from abcpy.statistics import Identity
statistics_calculator = Identity(degree = 1, cross = False)
# Redefine the statistics function
statistics_calculator.statistics = lambda x, f1=statistics_calculator.statistics: (f1(x)[:,indices_chosen])

XObserved = np.zeros(shape=(40,25))
Xtransformed = np.zeros(shape=(40,10))
for whichobs in range(40):
    XObserved[whichobs,:] = [np.hstack((np.array([0,20,60,120,300]).reshape(5,1),X[whichobs,:].reshape(4,5).transpose())).flatten().reshape(1,-1)][0]
    Xtransformed[whichobs, :] = statistics_calculator.statistics([XObserved[whichobs,:]])

# from abcpy.statistics import Identity
# statistics_calculator = Identity(degree = 1, cross = False)
# # Redefine the statistics function
# L = 1/np.std(Xtransformed, axis=0)
# statistics_calculator.statistics = lambda x, f1=statistics_calculator.statistics: (f1(x)[:,indices_chosen])*L
# for whichobs in range(40):
#     XObserved[whichobs,:] = [np.hstack((np.array([0,20,60,120,300]).reshape(5,1),X[whichobs,:].reshape(4,5).transpose())).flatten().reshape(1,-1)][0]
#     Xtransformed[whichobs, :] = statistics_calculator.statistics([XObserved[whichobs,:]])


# Consider 3 classes as volunteer, Patient Dialysé (non Diabètiqu + Diabètiqu) and  Patient BPCO (stable + Exacerbé)
# A, B, C = X[:7, :], X[7:23, :], X[23:, :]
# AI, BI, CI = np.ones(shape=(7,1)), np.ones(shape=(16,1)), np.ones(shape=(8,1))
# X = np.vstack((A, B, C))
# Y = np.vstack((AI, 2*BI, 3*CI)).reshape(-1,)
# n_clus = 3

n_clus = 3
A, B, C = Xtransformed[:16, :], Xtransformed[16:32, :], Xtransformed[32:, :]
AI, BI, CI = np.ones(shape=(16,1)), np.ones(shape=(16,1)), np.ones(shape=(8,1))
#A, B, C = Xtransformed[9:16, :], Xtransformed[16:32, :], Xtransformed[32:, :]
#AI, BI, CI = np.ones(shape=(7,1)), np.ones(shape=(16,1)), np.ones(shape=(8,1))
Y = np.vstack((AI, 2*BI, 3*CI)).reshape(-1,)
X = np.vstack((A, B, C))

# Hierarchical clustering using Euclidean
cluster = AgglomerativeClustering(n_clusters=n_clus, affinity='euclidean', linkage='ward')
cluster.fit_predict(X)
print('AR score of Euclidean distance: '+str(adjusted_rand_score(Y, cluster.labels_)))

# Hierachical clustering using Euclidean distance between the projections learned using Large Margin Nearest Neighbor (LMNN) algorithm
from metric_learn import LMNN, NCA, MLKR, MMC_Supervised
lmnn = LMNN(init='lda',k=6, min_iter=10000, max_iter=50000, convergence_tol=1e-6, learn_rate=1e-7, regularization=.5)
# lmnn.fit(X, Y)
# X_lmnn = lmnn.transform(X)
X_lmnn = lmnn.fit_transform(X, Y)
print(X_lmnn.shape)
cluster_lmnn = AgglomerativeClustering(n_clusters=n_clus, affinity='euclidean', linkage='ward')
cluster_lmnn.fit_predict(X_lmnn)
print(X_lmnn[0,:].shape)
# score[ind] = adjusted_rand_score(Y, cluster_lmnn.labels_)
# print(score[ind])
print('AR score of Euclidean distance btwn LMNN: '+str(adjusted_rand_score(Y, cluster_lmnn.labels_)))

colorarray = ['tomato', 'royalblue', 'limegreen']
fig, ax = plt.subplots(1)
ax.plot(X_lmnn[:16,0],X_lmnn[:16,1], '*', color=colorarray[0],  label='Healthy volunteers')
ax.plot(X_lmnn[16:32,0],X_lmnn[16:32,1], '*',  color=colorarray[1], label='Dialysis patients')
ax.plot(X_lmnn[32:,0],X_lmnn[32:,1], '*', color=colorarray[2], label='COPD patients')
ax.legend(fancybox=True, framealpha=0.5, fontsize=10)
ax.set_xlabel(r'$X_1$', fontsize=15)
ax.set_ylabel(r'$X_2$', fontsize=15)
ax.set_title('ProjectedSpace')
plt.tight_layout()
#for ind in range(10):
#    plt.plot(X_lmnn[ind,0], X_lmnn[ind,1], color = colorarray[int(Y[ind])])
plt.savefig('ProjectedSpace.pdf')

# nca = NCA(init='auto',max_iter=10000)
# nca.fit(X, Y)
# X_nca = nca.transform(X)
# cluster_nca = AgglomerativeClustering(n_clusters=n_clus, affinity='euclidean', linkage='ward')
# cluster_nca.fit_predict(X_nca)
# print('AR score of Euclidean distance btwn NCA: '+str(adjusted_rand_score(Y, cluster_nca.labels_)))
#
# mlkr = MLKR(init='identity')
# mlkr.fit(X, Y)
# X_mlkr = mlkr.transform(X)
# cluster_mlkr = AgglomerativeClustering(n_clusters=n_clus, affinity='euclidean', linkage='ward')
# cluster_mlkr.fit_predict(X_mlkr)
# print('AR score of Euclidean distance btwn MLKR: '+str(adjusted_rand_score(Y, cluster_mlkr.labels_)))
#
# mmc = MMC_Supervised(init='identity')
# mmc.fit(X, Y)
# X_mmc = mmc.transform(X)
# cluster_mmc = AgglomerativeClustering(n_clusters=n_clus, affinity='euclidean', linkage='ward')
# cluster_mmc.fit_predict(X_mmc)
# print('AR score of Euclidean distance btwn MMC: '+str(adjusted_rand_score(Y, cluster_mmc.labels_)))

