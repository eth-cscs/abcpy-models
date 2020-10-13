import pylab as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score
from abcpy.output import Journal
import numpy as np
import os.path
from scipy.stats import gaussian_kde

names = ['pAD','pAg','pT','pF','aT','v_z_AP','v_z_NAP']

margmax, meanpost, ylabel = [], [], []
for ind in range(40):
    whichobs = ind
    filename = 'apmcabc_obs_'+str(whichobs)+'.jrnl'
    if os.path.isfile(filename):
        journal = Journal.fromFile(filename)
        weights = np.concatenate(journal.get_weights())
        post1 = np.concatenate(journal.get_parameters()['pAD'])
        post2 = np.concatenate(journal.get_parameters()['pAg'])
        post3 = np.concatenate(journal.get_parameters()['pT'])
        post4 = np.concatenate(journal.get_parameters()['pF'])
        post5 = np.concatenate(journal.get_parameters()['aT'])
        post6 = np.concatenate(journal.get_parameters()['v_z_AP'])
        post7 = np.concatenate(journal.get_parameters()['v_z_NAP'])
        data = np.hstack((post1, post2, post3, post4, post5, post6, post7))
        meanpost.append(np.array([np.mean(np.concatenate(journal.get_parameters()[x])) for x in names]))
        #meanpost.append(np.array([journal.posterior_mean()[x] for x in names]))
        # marginal_max = np.zeros(shape=(len(names),))
        # for i, label in enumerate(names):
        #     xmin, xmax = min(data[:,i]), max(data[:,i])
        #     positions = np.linspace(xmin, xmax, 100)
        #     gaussian_kernel = gaussian_kde(data[:,i], weights=weights, bw_method=.4)
        #     AA = gaussian_kernel(positions)
        #     marginal_max[i] = positions[AA == max(AA)]
        # margmax.append(marginal_max)
        ylabel.append(Y[whichobs])

# #print(meanpost)
# meanpost = np.array(meanpost)
# #margmax = np.array(margmax)
# ylabel = np.array(ylabel)
#
# plt.figure()
# plt.plot(meanpost[10:16,0], meanpost[10:16,1], 'b*')
# plt.plot(meanpost[ylabel==2,0], meanpost[ylabel==2,1], 'r*')
# plt.savefig('scatter.eps')
# plt.show()
#
# print(meanpost[ylabel==1,:])
#
# contr_color1 = 'tomato'
# contr_color2 = 'red'
# tripl_color1 = 'royalblue'
# tripl_color2 = 'blue'
# sdml_color1 = 'limegreen'
# sdml_color2 = 'green'
# alpha = 0.5
#
# def set_colors(bp, color1, color2, alpha):
#     for box in bp['boxes']:
#         # change outline color
#         box.set(color=color1, alpha=alpha)
#         # change fill color
#         # box.set( facecolor = '#1b9e77' )
#
#     ## change color and linewidth of the whiskers
#     for whisker in bp['whiskers']:
#         whisker.set(color=color1, alpha=alpha)  # , linewidth=2)
#
#     ## change color and linewidth of the caps
#     for cap in bp['caps']:
#         cap.set(color=color1, linewidth=2, alpha=alpha)
#
#     ## change color and linewidth of the medians
#     for median in bp['medians']:
#         median.set(color=color2, linewidth=2)
#
#     ## change color and linewidth of the medians
#     for mean in bp['means']:
#         mean.set(color=color2, linewidth=2)
#
#     ## change the style of fliers and their fill
#     for flier in bp['fliers']:
#         flier.set(marker='.', alpha=alpha, markerfacecolor=color1, markeredgecolor=color1)
#
# namestitle = [r'$p_{Ad}$',r'$p_{Ag}$',r'$p_T$',r'$p_F$',r'$a_T$','v_z_AP','v_z_NAP']
# for k in range(7):
#     fig, ax = plt.subplots(1)
#     ax.set_title(namestitle[k], fontsize=20)
#     ax.boxplot([meanpost[ylabel==1,k],meanpost[ylabel==2,k]], labels = ['Healthy voluteers', 'Dialysis patients'])
#     plt.savefig('Boxplot_meanpost_'+str(k)+'.eps')
#
# namestitle = [r'$p_{Ad}$',r'$p_{Ag}$',r'$p_T$',r'$p_F$',r'$a_T$','v_z_AP','v_z_NAP']
# for k in range(7):
#     fig, ax = plt.subplots(1)
#     ax.set_title(namestitle[k], fontsize=20)
#     ax.boxplot([meanpost[ylabel==1,k],meanpost[ylabel==2,k], meanpost[ylabel==3,k]], labels = ['Healthy voluteers', 'Dialysis patients', 'Patients with COPD'])
#     plt.savefig('Boxplot_meanpost_'+str(k)+'.eps')
#
# namestitle = [r'$p_{Ad}$',r'$p_{Ag}$',r'$p_T$',r'$p_F$',r'$a_T$','v_z_AP','v_z_NAP']
# for k in range(7):
#     fig, ax = plt.subplots(1)
#     ax.set_title(namestitle[k], fontsize=20)
#     ax.boxplot([meanpost[ylabel==1,k],meanpost[ylabel==2,k], meanpost[ylabel==3,k]], labels = ['Healthy voluteers', 'Dialysis patients', 'Patients with COPD'])
#     plt.savefig('Boxplot_margmax_'+str(k)+'.eps')
#
# k = 6
# fig, ax = plt.subplots(1)
# bp_contr = ax.boxplot(meanpost[ylabel==1,k], notch=False, meanline=True, showfliers=True, patch_artist=True)
# set_colors(bp_contr, contr_color1, contr_color2, alpha)
# bp_tripl = ax.boxplot(meanpost[ylabel==2,k], notch=False, meanline=True, showfliers=True, patch_artist=True)
# set_colors(bp_tripl, tripl_color1, tripl_color2, alpha)
# bp_sdml = ax.boxplot(meanpost[ylabel==3,k], notch=False, meanline=True, showfliers=True, patch_artist=True)
# set_colors(bp_sdml, sdml_color1, sdml_color2, alpha)
# ax.legend([bp_contr["boxes"][0], bp_tripl["boxes"][0], bp_sdml["boxes"][0]],
#           ['Healthy voluteers', 'Dialysis patients', 'Patients with COPD'])
# plt.savefig('Boxplot_'+str(k)+'.eps')
# plt.show()
#
#
# # Hierarchical clustering using Euclidean
# cluster = AgglomerativeClustering(n_clusters=n_clus, affinity='euclidean', linkage='ward')
# cluster.fit_predict(meanpost[:,:4])
# print('AR score of Euclidean distance: '+str(adjusted_rand_score(ylabel, cluster.labels_)))
#
# colorarray = ['b', 'r', 'k']
# fig, ax = plt.subplots(1)
# for ind in range(23):
#     ax.plot(meanpost[ind,3],meanpost[ind,4], '*', color = colorarray[int(ylabel[ind]-1)])
# ax.legend(fancybox=True, framealpha=0.5, fontsize=10)
# ax.set_xlabel(r'$X_1$', fontsize=15)
# ax.set_ylabel(r'$X_2$', fontsize=15)
# plt.tight_layout()
# plt.show()
# for ind in range(10):
#    plt.plot(X_lmnn[ind,0], X_lmnn[ind,1], color = colorarray[int(Y[ind])])