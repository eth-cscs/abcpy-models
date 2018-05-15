import numpy as np
import scipy.sparse.linalg as linalg
import networkx as nx
import matplotlib.pyplot as plt

#==============================================================================
# Define NetSleuth Algorithm
#==============================================================================
def netsleuth(infectedi, G, beta):
    G = beta * G
    infected = infectedi
    if len(infected) <= 1:
        seedids = infected
    else:
        # take submatrix
        A = G[infected, :][:, infected]
        factor, factor_i = 1, 1
        DI = factor_i * np.diag(sum(A))
        TD = np.diag(sum(G)[infected])
        DNI = factor*(TD - DI)
        LA = DNI + DI - A

        # find smallest eig
        [l, u] = linalg.eigs(LA, k=1, which='SM',tol=0.01)
        u, l = np.abs(u), np.abs(l)
        u_copy = u
        ind = np.zeros(shape=(len(u),)).astype(int)
        for i in range(len(u)):
            val, ind[i] = np.amax(u_copy), int(np.argmax(u_copy))
            u_copy[ind[i]] = -99

        for i in range(len(ind)):
            ids = infected[int(ind[i])]
            if len(np.where(np.unique(infectedi) == ids)) != 0:
                seedids = ids
                break

    return seedids

#==============================================================================
# Different types of network (BA: Barabasi-Albert, ER: Erdos-Renyi, FB: Facebook Social Network,
# INRV: Indian Village contact Network) with node_no many nodes on the network. The infection_node
# is the true seed-node. (Choose one of the options)
#==============================================================================
#case, node_no, infection_start_point, const_figure_tick = 'ba', 100, 4, .69
case, node_no, infection_start_point, const_figure_tick = 'er', 100, 10, .39
#==============================================================================
# Load network
#==============================================================================
A = np.load('Networks/'+case+'_'+str(node_no)+'_network.npy')
network = nx.from_numpy_matrix(A)
# ==============================================================================
# Load simulated dataset
# ==============================================================================
y_obs = np.load('Results/SimpleContagion/' + case + '_' + str(node_no) + '_yobs_' + str(infection_start_point) + '.npy')
# ==============================================================================
# Compute seed node using NetSleuth
# ==============================================================================
beta = 0.3
seed_est = np.zeros(shape=(100,))
seed_est_pl = np.zeros(shape=(100,))
for ind in range(100):
    seed_ns = netsleuth(np.array(y_obs[ind][0][0]).astype(int), A, beta)
    print(str(ind)+': seed :'+str(seed_ns)+':shortest path length:'+str(nx.shortest_path_length(network,source=infection_start_point,target=seed_ns)))
    seed_est[ind] = seed_ns
    seed_est_pl[ind] = nx.shortest_path_length(network,source=infection_start_point,target=seed_ns)
# ==============================================================================
# Create Figure
# ==============================================================================
uniq_post_size = np.unique(seed_est_pl)
avg_post_prob = np.zeros(shape=(5,1))
for ind1 in range(len(uniq_post_size)):
    avg_post_prob[ind1] = np.sum([seed_est_pl==uniq_post_size[ind1]])/100
uniq_post_size = [1,2,3,4,5]
plt.figure()
bar_width, opacity, error_config = 0.35, 0.4, {'ecolor': '0.3'}
plt.bar(uniq_post_size, avg_post_prob, bar_width, alpha=opacity, color='gray',align='center',linewidth=0, edgecolor = 'None')
plt.xlabel(r'$\rho(\hat{n}_{\mathrm{SN}},n^0_{\mathrm{SN}})$', fontsize = 30)
plt.ylabel('Density', fontsize = 30)

x_ticks_val = [0,1,2,3,4,5]
y_ticks_val = np.round(np.arange(min(avg_post_prob),max(avg_post_prob)+0.1*(max(avg_post_prob)-min(avg_post_prob)),0.5*(max(avg_post_prob)-min(avg_post_prob))),2)
plt.yticks(y_ticks_val,fontsize=20)
plt.xticks(x_ticks_val,fontsize=20)
plt.xlim([-.5, 5.5])
plt.tight_layout()
plt.tick_params(top='off', right='off', which='both')
plt.savefig('Figure/'+case+'_SC_network_bias_sn_Netsleuth.eps', format='eps', dpi=1000)