import numpy as np
import networkx as nx

from ComplexContagion.Model import ComplexContagion
from ComplexContagion.Statistics import DiffusionIdentityStatistics
from ComplexContagion.Distance import SubsetDistance
from ComplexContagion.Kernel import DiffusionKernel
from ComplexContagion.Prior import DiffusionPrior
from ComplexContagion.Inference import SABCDiffusion

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
time_observed = np.arange(20, 120+1)
#==============================================================================
# Load network
#==============================================================================
A = np.load('Networks/'+case+'_'+str(node_no)+'_network.npy')
network = nx.from_numpy_matrix(A)
#==============================================================================
# Generate dataset
#==============================================================================
theta, gamma, seed = 0.7, 0.3, 1
# Complex Contagion process simulation
prior = DiffusionPrior(np.array([1,1,0]), network, network.nodes(),[0.0,0.0,0.0],[1.0,1.0,1.0],seed = 1)
for infection_start_point in [infection_node]:
    model = ComplexContagion(prior, network, time_observed, theta, gamma, [infection_start_point], seed)
    if case is 'ba' or case is 'er':
        y_obs = model.simulate(100)
    else:
        y_obs = model.simulate(1)
    np.save('Results/ComplexContagion/'+case+'_yobs_'+str(infection_start_point)+'.npy',y_obs)
#==============================================================================
# Inference
#==============================================================================
# Define Statistics and Distance
stat_calc = DiffusionIdentityStatistics(degree=1, cross=0)
dist_calc = SubsetDistance(stat_calc, network,sp_distance = 1)

for infection_start_point in [infection_node]:
     # Load simulated dataset
     y_obs = np.load('Results/ComplexContagion/'+case+'_yobs_'+str(infection_start_point)+'.npy')
     print('Infection Start point:' + str(infection_start_point))
     for ind in range(len(y_obs)):
         kernel = DiffusionKernel(network, mean = 0, cov = 1, seed_node = 1, seed = 1)
         #Redefine prior given observed dataset       
         prior = DiffusionPrior(np.array([1,1,0]), network, np.array(y_obs[ind][0][0]),[0.0,0.0,0.0],[1.0,1.0,1.0],seed = 1)
         #Initiate model with new prior
         model = ComplexContagion(prior, network, time_observed, theta = .3, gamma = .5, infection_start_point = [0], seed = 1)
         # Sampling using SABC       
         sampler_sabc = SABCDiffusion(model, dist_calc, kernel, backend, seed = 1)
         step, epsilon_sabc, n_samples, n_samples_per_param, ar_cutoff, beta, delta, v = 200, np.array([40.0]), 1000, 1, 0.0001, 2, 0.2, 0.3
         journal_sabc = sampler_sabc.sample([y_obs[ind]], step, epsilon_sabc, n_samples, n_samples_per_param, beta, delta, v, ar_cutoff, resample=None, n_update=None, adaptcov=1, full_output=1)
         journal_sabc.save('Results/ComplexContagion/'+case+'_'+str(node_no)+'_joint_SABC_'+str(infection_start_point)+'_'+str(ind)+'.jrnl')
        
     del A, network, y_obs, kernel, prior, model, sampler_sabc, journal_sabc
