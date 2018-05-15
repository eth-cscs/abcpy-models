import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from abcpy.output import Journal
from scipy.stats import gaussian_kde

#==============================================================================
# Different types of network (BA: Barabasi-Albert, ER: Erdos-Renyi, FB: Facebook Social Network,
# INRV: Indian Village contact Network) with node_no many nodes on the network. The infection_node
# is the true seed-node. (Choose one of the options)
#==============================================================================
#case, node_no, infection_start_point, ind = 'ba', 100, 4, 21
#case, node_no, infection_start_point, ind = 'er', 100, 10, 9
#case, node_no, infection_start_point, ind = 'inrv', 354, 70, 0
case, node_no, infection_start_point, ind = 'fb', 4039, 2000, 0
#==============================================================================
# Load network
#==============================================================================
A = np.load('Networks/'+case+'_'+str(node_no)+'_network.npy')
network = nx.from_numpy_matrix(A)
#==============================================================================
# Load the inferred posterior distribution
#==============================================================================
jj = Journal.fromFile('Results/ComplexContagion/'+case+'_'+str(node_no)+'_joint_SABC_'+str(infection_start_point)+'_'+str(ind)+'.jrnl')
samples = jj.get_parameters()
#==============================================================================
# Compute Bayes Estimate
#==============================================================================
from BayesEstimate import BayesEstimateComplexContagion
BE = BayesEstimateComplexContagion(network, samples, 10)
Bayes_estimate = BE.minimize()
print(Bayes_estimate)
print(nx.shortest_path_length(network, infection_start_point, Bayes_estimate[2]))
#==============================================================================
# Compute and plot the posterior distribution for the dataset at 'ind'
#==============================================================================
xmin, xmax = min(samples[:,0]), max(samples[:,0])
ymin, ymax = min(samples[:,1]), max(samples[:,1])
X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([samples[:,0], samples[:,1]])
kernel = gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)
plt.figure()
plt.plot(samples[:,0],samples[:,1],'.',color = '0.8',alpha = 0.2)
CS = plt.contour(X, Y, Z, 6, colors = 'k')
CS = plt.contour(X, Y, Z, 6, colors = 'k',linestyles='solid')
plt.clabel(CS, fontsize=9, inline=1)
plt.plot(.7, .3, 'k',linestyle='solid', markersize=4,label = r'$P(\theta|x_0)$')
plt.plot(0.7, 0.3, '+k',label= r'$\theta^0$', markersize = 30)
plt.plot(Bayes_estimate[0], Bayes_estimate[1], 'x',color='k',label= r'$\hat{\theta}$', markersize = 30)
plt.plot([0.7,0.7],[0,1],'k',linestyle='solid')
plt.plot([0,1],[0.3,0.3],'k',linestyle='solid')
plt.xlabel(r'$\beta$', fontsize=30)
plt.ylabel(r'$\gamma$', fontsize=30)
plt.xticks(np.arange(0,1.1,.5),fontsize=20)
plt.yticks(np.arange(0,1.1,.5),fontsize=20)
plt.xlim([0,1])
plt.ylim([0,1])
plt.tight_layout()
plt.legend(loc='upper left', frameon=False, numpoints=1,fontsize=25)
plt.savefig('Figure/'+case+'_CC_posterior_distribution_'+str(ind)+'.eps', format='eps', dpi=1000)


A = np.unique(samples[:,2])
H = nx.subgraph(network,A)
post_tmp = np.zeros(shape=(len(A),))
post_size = np.ones(shape=(len(A),))
for ind1 in range(len(H.nodes())):
    post_tmp[ind1] = sum(samples[:, 2] == H.nodes()[ind1]) / samples.shape[0]
    post_size[ind1] = nx.shortest_path_length(network, H.nodes()[ind1], infection_start_point)
uniq_post_size = np.unique(post_size)
sum_post_prob = np.zeros(shape=(6,1))
max_post_prob = np.zeros(shape=(6,1))
for ind1 in range(len(uniq_post_size)):
    sum_post_prob[ind1] = np.sum(post_tmp[post_size == uniq_post_size[ind1]])
    max_post_prob[ind1] = np.max(post_tmp[post_size == uniq_post_size[ind1]])
uniq_post_size = [0,1,2,3,4,5]

plt.figure()
bar_width, opacity, error_config = 0.35, 0.4, {'ecolor': '0.3'}
plt.bar(uniq_post_size, sum_post_prob, bar_width, alpha=1, color='gray',align='center',linewidth=0, edgecolor = 'None')
plt.xlabel(r'$\rho(n_{\mathrm{SN}},n^0_{\mathrm{SN}})$', fontsize = 30)
plt.ylabel('Density', fontsize = 30)
y_ticks_val = np.round(np.arange(min(sum_post_prob),max(sum_post_prob)+0.1*(max(sum_post_prob)-min(sum_post_prob)),0.5*(max(sum_post_prob)-min(sum_post_prob))),2)
x_ticks_val = np.round(np.arange(min(uniq_post_size),max(uniq_post_size)+0.1*(max(uniq_post_size)-min(uniq_post_size)),0.2*(max(uniq_post_size)-min(uniq_post_size))),2)
plt.yticks(y_ticks_val,fontsize=20)
plt.xticks(x_ticks_val,fontsize=20)
plt.tight_layout()
plt.savefig('Figure/'+case+'_CC_posterior_distribution_'+str(ind)+'_'+'sn.eps', format='eps', dpi=1000)

plt.figure()
bar_width, opacity, error_config = 0.35, 0.4, {'ecolor': '0.3'}
plt.bar(uniq_post_size, max_post_prob, bar_width, alpha=1, color='gray',align='center',linewidth=0, edgecolor = 'None')
plt.xlabel(r'$\rho(n_{\mathrm{SN}},n^0_{\mathrm{SN}})$', fontsize = 30)
plt.ylabel('Density', fontsize = 30)
y_ticks_val = np.round(np.arange(min(max_post_prob),max(max_post_prob)+0.1*(max(max_post_prob)-min(max_post_prob)),0.5*(max(max_post_prob)-min(max_post_prob))),2)
x_ticks_val = np.round(np.arange(min(uniq_post_size),max(uniq_post_size)+0.1*(max(uniq_post_size)-min(uniq_post_size)),0.2*(max(uniq_post_size)-min(uniq_post_size))),2)
plt.yticks(y_ticks_val,fontsize=20)
plt.xticks(x_ticks_val,fontsize=20)
plt.tight_layout()
plt.savefig('Figure/'+case+'_CC_posterior_distribution_'+str(ind)+'_'+'sn_max.eps', format='eps', dpi=1000)