import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from abcpy.output import Journal
from scipy.stats import gaussian_kde

#==============================================================================
# Different types of network (BA: Barabasi-Albert, ER: Erdos-Renyi, FB: Facebook Social Network,
# INRV: Indian Village contact Network) with node_no many nodes on the network. The infection_node
# is the true seed-node.
#==============================================================================
case, node_no, infection_start_point, ind = 'ba', 100, 4, 27
#case, node_no, infection_start_point, ind = 'er', 100, 10, 17
#case, node_no, infection_start_point, ind = 'inrv', 354, 70, 0
#case, node_no, infection_start_point, ind = 'fb', 4039, 2000, 0
#==============================================================================
# Load network
#==============================================================================
A = np.load('Networks/'+case+'_'+str(node_no)+'_network.npy')
network = nx.from_numpy_matrix(A)
#==============================================================================
# Load the inferred posterior distribution
#==============================================================================
jj = Journal.fromFile('Results/SimpleContagion/'+case+'_'+str(node_no)+'_joint_SABC_'+str(infection_start_point)+'_'+str(ind)+'.jrnl')
samples = jj.get_parameters()
#==============================================================================
# Compute Bayes Estimate
#==============================================================================
from BayesEstimate import BayesEstimateSimpleContagion
BE = BayesEstimateSimpleContagion(network, samples)
Bayes_estimate = BE.compute()
print(Bayes_estimate)
print(nx.shortest_path_length(network, infection_start_point, Bayes_estimate[1]))
#==============================================================================
# Compute and plot the posterior distribution for the dataset at 'ind'
#==============================================================================
theta = 0.3
xmin, xmax = 0, 1
positions = np.linspace(xmin, xmax, samples.shape[0])
gaussian_kernel = gaussian_kde(samples[:,0].reshape(samples.shape[0],),bw_method = .8)
values = gaussian_kernel(positions)
plt.figure()
plt.plot(positions,gaussian_kernel(positions), label = r'$P(\theta|x_0)$',color='k',linestyle='solid')
plt.plot([theta, theta],[0.0, np.max(values)+0.5], label = r'$\theta^0$',color ='k',linestyle='dashed')
plt.plot([np.mean(samples[:,0]), np.mean(samples[:,0])],[0.0, np.max(values)+0.5], label = r'$\hat{\theta}$',color ='k',linestyle='dashdot')
#plt.plot([Bayes_estimate_p[infection_node,0], Bayes_estimate_p[infection_node,0]],[0.0, np.max(values)+0.5], label = r'$\hat{\theta}$',color ='k',linestyle='dashdot')
plt.xlabel(r'$\theta$', fontsize = 30)
plt.ylabel('Density', fontsize = 30)
plt.ylim([0.0, np.max(values)+0.5])
plt.xlim([0,1])
plt.xticks(np.arange(0,1.1,.5),fontsize=20)
plt.yticks(np.round(np.arange(min(values),max(values)+1,0.5*(max(values)-min(values))),1),fontsize=20)
plt.legend(loc='best', frameon=False, numpoints=1, fontsize=25)
plt.tight_layout()
plt.savefig('Figure/'+case+'_SC_posterior_distribution_'+str(ind)+'_theta.eps', format='eps', dpi=1000)
plt.close()

A = np.unique(samples[:,1])
H = nx.subgraph(network,A)
post_tmp = np.zeros(shape=(len(A),))
post_size = np.ones(shape=(len(A),))
for ind1 in range(len(H.nodes())):
    post_tmp[ind1] = sum(samples[:, 1] == H.nodes()[ind1]) / samples.shape[0]
    post_size[ind1] = nx.shortest_path_length(network, H.nodes()[ind1], infection_start_point)
uniq_post_size = np.unique(post_size)
sum_post_prob = np.zeros(shape=(6,1))
max_post_prob = np.zeros(shape=(6,1))
for ind1 in range(len(uniq_post_size)):
    sum_post_prob[ind1] = np.sum(post_tmp[post_size==uniq_post_size[ind1]])
    max_post_prob[ind1] = np.max(post_tmp[post_size == uniq_post_size[ind1]])
uniq_post_size = [0,1,2,3,4,5]

plt.figure()
bar_width, opacity, error_config = 0.35, 0.4, {'ecolor': '0.3'}
plt.bar(uniq_post_size, sum_post_prob, bar_width, alpha=opacity, color='gray',align='center',linewidth=0, edgecolor = 'None')
plt.xlabel(r'$\rho(n_{\mathrm{SN}},n^0_{\mathrm{SN}})$', fontsize = 30)
plt.ylabel('Density', fontsize = 30)
y_ticks_val = np.round(np.arange(min(sum_post_prob),max(sum_post_prob)+0.1*(max(sum_post_prob)-min(sum_post_prob)),0.5*(max(sum_post_prob)-min(sum_post_prob))),2)
x_ticks_val = np.round(np.arange(min(uniq_post_size),max(uniq_post_size)+0.1*(max(uniq_post_size)-min(uniq_post_size)),0.2*(max(uniq_post_size)-min(uniq_post_size))),2)
plt.yticks(y_ticks_val,fontsize=20)
plt.xticks(x_ticks_val,fontsize=20)
plt.tight_layout()
plt.savefig('Figure/'+case+'_SC_posterior_distribution_'+str(ind)+'_sn.eps', format='eps', dpi=1000)
plt.close()

plt.figure()
bar_width, opacity, error_config = 0.35, 0.4, {'ecolor': '0.3'}
plt.bar(uniq_post_size, max_post_prob, bar_width, alpha=opacity, color='gray',align='center',linewidth=0, edgecolor = 'None')
plt.xlabel(r'$\rho(n_{\mathrm{SN}},n^0_{\mathrm{SN}})$', fontsize = 30)
plt.ylabel('Density', fontsize = 30)
y_ticks_val = np.round(np.arange(min(max_post_prob),max(max_post_prob)+0.1*(max(max_post_prob)-min(max_post_prob)),0.5*(max(max_post_prob)-min(max_post_prob))),2)
x_ticks_val = np.round(np.arange(min(uniq_post_size),max(uniq_post_size)+0.1*(max(uniq_post_size)-min(uniq_post_size)),0.2*(max(uniq_post_size)-min(uniq_post_size))),2)
plt.yticks(y_ticks_val,fontsize=20)
plt.xticks(x_ticks_val,fontsize=20)
plt.tight_layout()
plt.savefig('Figure/'+case+'_SC_posterior_distribution_'+str(ind)+'_sn_max.eps', format='eps', dpi=1000)
plt.close()