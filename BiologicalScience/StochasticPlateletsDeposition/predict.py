# import logging
# logging.basicConfig(level=logging.DEBUG)
import matplotlib.pyplot as plt
import numpy as np
from abcpy.continuousmodels import Uniform
from abcpy.output import Journal
from Model import PlateletDeposition
from abcpy.statisticslearning import TripletDistanceLearning, SemiautomaticNN
from statistic import IdentityChosen
from abcpy.statistics import NeuralEmbedding
# ########### Read Observed data ################
from numpy import genfromtxt
my_data = genfromtxt('AllResults.csv', delimiter=',')
###############################################
simulate = False
plot = True
analyse = False

# Define backend
from abcpy.backends import BackendMPI, BackendDummy
backend = BackendDummy()
#backend = BackendMPI()
########### Define Graphical Model ############
########### Define Graphical Model ############
# True parameter 89.0, 76.0, 2.49, 7e-3, 7.7, 6e-3, 8e-4
XObserved = np.array([0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.11600000e+05, 5.96585366e+03,
                      2.00000000e+01, 3.90572997e+03, 1.59328549e+01, 1.38943902e+05, 0.00000000e+00,
                      6.00000000e+01, 3.42727305e+03, 2.80052570e+01, 8.57585366e+04, 0.00000000e+00,
                      1.20000000e+02, 2.33523014e+03, 7.57715388e+01, 4.25231707e+04, 0.00000000e+00,
                      3.00000000e+02, 1.74166329e+02, 2.46413793e+03, 5.15975610e+03, 0.00000000e+00])
obsdata = [np.array(XObserved).reshape(1, -1)]
noAP, noNAP, SR_x = int(my_data[0, 16]), int(my_data[0, 11]), float(my_data[0, 21])
# Define the indices used as features to remove timestamp etc. and values with NA
AllIndices = list(np.arange(0, 25, 1))
NANTimeIndices = [0, 5, 10, 15, 20] + list(np.arange(0, 5, 1)) + [4, 9, 14, 19, 24]
InformativeIndices = [item for item in AllIndices if item not in NANTimeIndices]
## Summary statistics taking informative indices
identity = IdentityChosen(InformativeIndices=InformativeIndices, degree=1,
                          cross=False)  # to apply before computing the statistics
## Random values
# The parameters considered random and we want to infer
pAd = Uniform([[5], [150]], name='pAD')
pAg = Uniform([[5], [150]], name='pAg')
pT = Uniform([[0.1], [10.0]], name='pT')
pF = Uniform([[0.1e-3], [9.0e-3]], name='pF')
aT = Uniform([[0], [10]], name='aT')
v_z_AP = Uniform([[1.0e-3], [9.0e-3]], name='v_z_AP')
v_z_NAP = Uniform([[1.0e-4], [9.0e-4]], name='v_z_NAP')
PD = PlateletDeposition([noAP, noNAP, SR_x, pAd, pAg, pT, pF, aT, v_z_AP, v_z_NAP], name='PD')

if simulate:
    from Model import DrawFromPosterior

    dp = DrawFromPosterior([PD], backend)
    print('Prediction')
    parameters, simulations = dp.sample(journal_file='Journals/sabc_obs_fake_multiply.jrnl')
    np.savez('Journals/fake_multiply_predict.npz', parameters=parameters, simulations=simulations)

    parameters, simulations = dp.sample(journal_file='Journals/sabc_obs_fake_seminn.jrnl')
    np.savez('Journals/fake_seminn_predict.npz', parameters=parameters, simulations=simulations)

    parameters, simulations = dp.sample(journal_file='Journals/sabc_obs_fake_triplet.jrnl')
    np.savez('Journals/fake_triplet_predict.npz', parameters=parameters, simulations=simulations)

if plot:
    multiply_predict = np.load('Journals/fake_multiply_predict.npz')
    seminn_predict = np.load('Journals/fake_seminn_predict.npz')
    triplet_predict = np.load('Journals/fake_triplet_predict.npz')
    simulations_multiply = multiply_predict['simulations']
    simulations_seminn = seminn_predict['simulations']
    simulations_triplet = triplet_predict['simulations']
    x_val = [20, 60, 120, 300]
    InformativeIndices = [[6, 11, 16, 21], [7, 12, 17, 22], [8, 13, 18, 23]]
    nameslabel = [r'$\mathbb{N}_{agg-clust} $', r'$\mathcal{S}_{agg-clust}}$', r'$\mathbb{N}_{platelet}$']
    filename = ['naggclust', 'saggclust', 'np']
    fig, axs = plt.subplots(1, 3, figsize=(30,10))  # 1
    for ind, ax in enumerate(axs.ravel()):  # 2
        #ax.set_title("Plot #{}".format(ind))
        ax.set_title(nameslabel[ind], fontsize=35)
        ax.plot(x_val, obsdata[0][0, InformativeIndices[ind]], color="black", marker=".", label=r'$\mathcal{x}^0$')

        ax.plot(x_val, np.quantile(simulations_seminn[:, InformativeIndices[ind]], .025, axis=0),color="r",linestyle="dotted", label='SASL')
        ax.plot(x_val, np.quantile(simulations_seminn[:, InformativeIndices[ind]], .975, axis=0),color="r",linestyle="dotted")
        ax.fill_between(x_val, np.quantile(simulations_seminn[:, InformativeIndices[ind]], .025, axis=0),
                        np.quantile(simulations_seminn[:, InformativeIndices[ind]], 0.975, axis=0), color='red', alpha=.1)
        ax.plot(x_val, np.quantile(simulations_seminn[:, InformativeIndices[ind]], .5, axis=0), color='red',marker=".",linestyle="dashed")

        ax.plot(x_val, np.quantile(simulations_multiply[:, InformativeIndices[ind]], .025, axis=0),color="b",linestyle="dotted", label='DSSL')
        ax.plot(x_val, np.quantile(simulations_multiply[:, InformativeIndices[ind]], .975, axis=0),color="b",linestyle="dotted")
        ax.fill_between(x_val, np.quantile(simulations_multiply[:, InformativeIndices[ind]], .025, axis=0),
                        np.quantile(simulations_multiply[:, InformativeIndices[ind]], 0.975, axis=0), color='b', alpha=.1)
        ax.plot(x_val, np.quantile(simulations_multiply[:, InformativeIndices[ind]], .5, axis=0), color='b',marker=".",linestyle="dashed")

        ax.plot(x_val, np.quantile(simulations_triplet[:, InformativeIndices[ind]], .025, axis=0),color="g",linestyle="dotted", label='TLSL')
        ax.plot(x_val, np.quantile(simulations_triplet[:, InformativeIndices[ind]], .975, axis=0),color="g",linestyle="dotted")
        ax.fill_between(x_val, np.quantile(simulations_triplet[:, InformativeIndices[ind]], .025, axis=0),
                        np.quantile(simulations_triplet[:, InformativeIndices[ind]], 0.975, axis=0), color='g', alpha=.1)
        ax.plot(x_val, np.quantile(simulations_triplet[:, InformativeIndices[ind]], .5, axis=0), color='g',marker=".",linestyle="dashed")
        leg = ax.legend(fontsize=30, fancybox=True)
        leg.get_frame().set_alpha(0)
        ax.set_xlabel("Time " + r"$(s)$", fontsize=30)
        ax.set_xticks(x_val)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax.tick_params(axis='both', labelsize=30)

        plt.tight_layout()
    plt.savefig('FinalFigures/Fig7.pdf')

if analyse:
    multiply_predict = np.load('Journals/fake_multiply_predict.npz')
    seminn_predict = np.load('Journals/fake_seminn_predict.npz')
    triplet_predict = np.load('Journals/fake_triplet_predict.npz')
    simulations_multiply = multiply_predict['simulations']
    simulations_seminn = seminn_predict['simulations']
    simulations_triplet = triplet_predict['simulations']

    InformativeIndices = [[6, 11, 16, 21], [7, 12, 17, 22], [8, 13, 18, 23]]
    nameslabel = [r'$\mathbb{N}_{agg-clust} $', r'$\mathcal{S}_{agg-clust}}$', r'$\mathbb{N}_{platelet}$']

    def score(data, dataobs, indices=None):
        print(data.shape, dataobs.shape)
        if indices is None:
            indices = [6, 11, 16, 21] + [7, 12, 17, 22] + [8, 13, 18, 23]
        score_1, score_2 = [], []
        for ind1 in range(data.shape[0]):
            score_1.append(2 * np.linalg.norm(dataobs[0,indices] - data[ind1,indices]))
        for ind1 in range(data.shape[0]):
            for ind2 in range(ind1,data.shape[0]):
                score_2.append(- np.linalg.norm(data[ind1,indices] - data[ind2, indices]))
        score = np.mean(score_1) + np.mean(score_2)
        return score

    print('All time series')
    print('score for multiply: '+str(score(simulations_multiply, obsdata[0], indices=None)))
    print('score for seminn: '+str(score(simulations_seminn, obsdata[0], indices=None)))
    print('score for triplet: '+str(score(simulations_triplet, obsdata[0], indices=None)))

    print('time series: '+nameslabel[0])
    print('score for multiply: '+str(score(simulations_multiply, obsdata[0], indices=[6, 11, 16, 21])))
    print('score for seminn: '+str(score(simulations_seminn, obsdata[0], indices=[6, 11, 16, 21])))
    print('score for triplet: '+str(score(simulations_triplet, obsdata[0], indices=[6, 11, 16, 21])))

    print('time series: '+nameslabel[1])
    print('score for multiply: '+str(score(simulations_multiply, obsdata[0], indices=[7, 12, 17, 22])))
    print('score for seminn: '+str(score(simulations_seminn, obsdata[0], indices=[7, 12, 17, 22])))
    print('score for triplet: '+str(score(simulations_triplet, obsdata[0], indices=[7, 12, 17, 22])))

    print('time series: '+nameslabel[0])
    print('score for multiply: '+str(score(simulations_multiply, obsdata[0], indices=[8, 13, 18, 23])))
    print('score for seminn: '+str(score(simulations_seminn, obsdata[0], indices=[8, 13, 18, 23])))
    print('score for triplet: '+str(score(simulations_triplet, obsdata[0], indices=[8, 13, 18, 23])))