import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from abcpy.output import Journal
from scipy.stats import gaussian_kde
from BayesEstimate import BayesEstimateComplexContagion

#==============================================================================
# Different types of network (BA: Barabasi-Albert, ER: Erdos-Renyi, FB: Facebook Social Network,
# INRV: Indian Village contact Network) with node_no many nodes on the network. The infection_node
# is the true seed-node.
#==============================================================================
#case, node_no, infection_start_point = 'ba', 100, 4
case, node_no, infection_start_point = 'er', 100, 10
#==============================================================================
# Load network
#==============================================================================
A = np.load('Networks/'+case+'_'+str(node_no)+'_network.npy')
network = nx.from_numpy_matrix(A)
# ==============================================================================
# Computation of Bayes Estimates
# ==============================================================================
Bayes_estimate_sn = np.zeros(shape=(100,))
Bayes_estimate_theta = np.zeros(shape=(100,2))
for ind in np.arange(100):
    print(ind)
    # ==============================================================================
    # Load the inferred posterior distribution
    # ==============================================================================
    jj = Journal.fromFile('Results/ComplexContagion/' + case + '_' + str(node_no) + '_joint_SABC_' + str(
        infection_start_point) + '_' + str(ind) + '.jrnl')
    samples = jj.get_parameters()
    # ==============================================================================
    # Compute Bayes Estimate
    # ==============================================================================
    BE = BayesEstimateComplexContagion(network, samples,10)
    Bayes_estimate = BE.minimize()
    print(Bayes_estimate)
    Bayes_estimate_theta[ind,:] = Bayes_estimate[:1]
    Bayes_estimate_sn[ind] = Bayes_estimate[2]
# ==============================================================================
# Plot the density of Bayes Estimates
# ==============================================================================
xmin, xmax = 0, 1
ymin, ymax = 0, 1
X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([Bayes_estimate_theta[:,0], Bayes_estimate_theta[:,1]])
kernel = gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)
plt.figure()
plt.plot(Bayes_estimate_theta[:,0],Bayes_estimate_theta[:,1],'.',color = '0.8',alpha=0.2)
CS = plt.contour(X, Y, Z, 6, colors = 'k',linestyles='solid')
plt.clabel(CS, fontsize=9, inline=1)
plt.plot(0.7, 0.3, '+k',label= r'$\theta^0$', markersize = 15)
plt.plot([0.7,0.7],[0,1],'k', linestyle='dashed')
plt.plot([0,1],[0.3,0.3],'k', linestyle='dashed')
plt.xlabel(r'$\hat{\beta}$', fontsize=30)
plt.ylabel(r'$\hat{\gamma}$', fontsize=30)
plt.xticks(np.arange(0,1.1,.5),fontsize=20)
plt.yticks(np.arange(0,1.1,.5),fontsize=20)
plt.legend(loc='best', frameon=False, numpoints=1, fontsize=25)
plt.tight_layout()
plt.savefig('Figure/'+case+'_CC_network_bias_beta_gamma.eps', format='eps', dpi=1000)

A = np.unique(Bayes_estimate_sn)
H = nx.subgraph(network,A)
post_tmp = np.zeros(shape=(len(A),))
post_size = np.ones(shape=(len(A),))
for ind1 in range(len(H.nodes())):
    post_tmp[ind1] = sum(Bayes_estimate_sn == H.nodes()[ind1]) / 100
    post_size[ind1] = nx.shortest_path_length(network, H.nodes()[ind1], infection_start_point)
uniq_post_size = np.unique(post_size)
avg_post_prob = np.zeros(shape=(6,1))
for ind1 in range(len(uniq_post_size)):
    avg_post_prob[ind1] = np.sum(post_tmp[post_size==uniq_post_size[ind1]])
uniq_post_size = [0,1,2,3,4,5]
plt.figure()
bar_width, opacity, error_config = 0.35, 0.4, {'ecolor': '0.3'}
plt.bar(uniq_post_size, avg_post_prob, bar_width, alpha=opacity, color='gray',align='center',linewidth=0, edgecolor = 'None')
plt.xlabel(r'$\rho(\hat{n}_{\mathrm{SN}},n^0_{\mathrm{SN}})$', fontsize = 30)
plt.ylabel('Density', fontsize = 30)
y_ticks_val = np.round(np.arange(min(avg_post_prob),max(avg_post_prob)+0.1*(max(avg_post_prob)-min(avg_post_prob)),0.5*(max(avg_post_prob)-min(avg_post_prob))),2)
x_ticks_val = np.round(np.arange(min(uniq_post_size),max(uniq_post_size)+0.1*(max(uniq_post_size)-min(uniq_post_size)),0.2*(max(uniq_post_size)-min(uniq_post_size))),2)
plt.yticks(y_ticks_val,fontsize=20)
plt.xticks(x_ticks_val,fontsize=20)
plt.tight_layout()
plt.savefig('Figure/'+case+'_CC_network_bias_sn.eps', format='eps', dpi=1000)