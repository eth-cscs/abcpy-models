import numpy as np
import networkx as nx

from SimpleContagion.Model import SimpleContagion
from SimpleContagion.Statistics import DiffusionIdentityStatistics
from SimpleContagion.Distance import SubsetDistance
from SimpleContagion.Kernel import DiffusionKernel
from SimpleContagion.Prior import DiffusionPrior
from SimpleContagion.Inference import SABCDiffusion


#==============================================================================
# Choose the appropriate Backend for Parallelization
from abcpy.backends import BackendDummy
backend = BackendDummy()
#from abcpy.backends import BackendMPI as Backend
#backend = Backend()
#==============================================================================
# Different types of network (BA: Barabasi-Albert, ER: Erdos-Renyi, FB: Facebook Social Network,
# INRV: Indian Village contact Network) with node_no many nodes on the network. The infection_node
# is the true seed-node. (Choose one of the options)
#==============================================================================
case, node_no, infection_node = 'ba', 100, 4
#case, node_no, infection_node = 'er', 100, 10
#case, node_no, infection_node = 'inrv', 354, 70
#case, node_no, infection_node = 'fb', 4039, 2000
#==============================================================================
# Time observed
time_observed = np.arange(20, 70+1)
#==============================================================================
# Load network
#==============================================================================
A = np.load('Networks/'+case+'_'+str(node_no)+'_network.npy')
network = nx.from_numpy_matrix(A)
#==============================================================================
# Generate dataset
#==============================================================================
theta, seed = 0.3, 1
prior = DiffusionPrior(np.array([1,0]), network, network.nodes(),[0.0,0.0],[1.0,1.0],seed = 1)
for infection_start_point in [infection_node]:   
    model = SimpleContagion(prior, network, time_observed, theta, infection_start_point, seed)
    if case is 'ba' or case is 'er':
        y_obs = model.simulate(100)
    else:
        y_obs = model.simulate(1)
    np.save('Results/SimpleContagion/'+case+'_'+str(node_no)+'_yobs_'+str(infection_start_point)+'.npy',y_obs)
#==============================================================================
# Inference
#==============================================================================
# Define Statistics and Distance
stat_calc = DiffusionIdentityStatistics(degree=1, cross=0)
dist_calc = SubsetDistance(stat_calc, network)

for infection_start_point in [infection_node]:
    # Load simulated dataset
    y_obs = np.load('Results/SimpleContagion/'+case+'_'+str(node_no)+'_yobs_'+str(infection_start_point)+'.npy')
    print('Infection Start point:' + str(infection_start_point))
    for ind in range(len(y_obs)):
        print(ind)           
        kernel = DiffusionKernel(network, mean = 0, var = 1, seed_node = 1, seed = 1)
        #Redefine prior given observed dataset       
        prior = DiffusionPrior(np.array([1,0]), network, np.array(y_obs[ind][0][0]),[0.0,0.0],[1.0,1.0],seed = 1)
        #Initiate model with new prior
        model = SimpleContagion(prior, network, time_observed, theta, infection_start_point, seed = None)
        # Sampling using SABC       
        sampler_sabc = SABCDiffusion(model, dist_calc, kernel, backend, seed = 1)
        step, epsilon_sabc, n_samples, n_samples_per_param, ar_cutoff, beta, delta, v = 200, np.array([40.0]), 1000, 1, 0.0001, 2, 0.2, 0.3
        journal_sabc = sampler_sabc.sample([y_obs[ind]], step, epsilon_sabc, n_samples, n_samples_per_param, beta, delta, v, ar_cutoff, resample=None, n_update=None, adaptcov=1, full_output=1)
        journal_sabc.save('Results/SimpleContagion/'+case+'_'+str(node_no)+'_joint_SABC_'+str(infection_start_point)+'_'+str(ind)+'.jrnl')
    
    del A, network, y_obs, kernel, prior, model, sampler_sabc, journal_sabc
