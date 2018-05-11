import numpy as np
import networkx as nx

from ComplexContagion.Model import ComplexContagion
from ComplexContagion.Statistics import DiffusionIdentityStatistics
from ComlexContagion.Distance import SubsetDistance
from ComplexContagion.Kernel import DiffusionKernel
from ComplexContagion.Prior import DiffusionPrior
from ComplexContagion.Inference import SABCDiffusion


from Diffusion.Diffusion_SC import Network_generate
import time


from abcpy.backends import BackendDummy
backend = BackendDummy()

#from abcpy.backends import BackendMPI as Backend
#backend = Backend()

# Different cases
case, node_no, infection_node = 'fb', 4039, 2000 # other cases 'ba', 'inrv', 'fb', 'tw', 'gog'
# Time observed
time_observed = np.arange(20, 120+1)
# Load network
A = np.load('Diffusion/Results/Contagion/'+case+'_'+str(node_no)+'_network.npy')
network = nx.from_numpy_matrix(A)

# Generate Network
#==============================================================================
#generator = Network_generate(node_no,1)
#network = generator.erdos_renyi(0.05)
#network.nodes()
#A = nx.to_numpy_matrix(network)
#np.save('Diffusion/Results/Contagion/'+case+'_'+str(node_no)+'_network.npy',A)
# 
#==============================================================================

# Generate dataset
#==============================================================================
theta, gamma, seed = 0.7, 0.3, 1
# Complex Contagion process simulation
prior = DiffusionPrior(np.array([1,1,0]), network, network.nodes(),[0.0,0.0,0.0],[1.0,1.0,1.0],seed = 1)
for infection_start_point in [infection_node]:#range(node_no):   
    model = ComplexContagion(prior, network, time_observed, theta, gamma, [infection_start_point], seed)
    y_obs = model.simulate(1)       
    np.save('Diffusion/Results/Contagion/'+case+'_yobs_'+str(infection_start_point)+'.npy',y_obs)
# Inference
#==============================================================================
# 
# Import Statistics and Distance
#==============================================================================
stat_calc = DiffusionIdentityStatistics(degree=1, cross=0)
dist_calc = SubsetDistance(stat_calc, network,sp_distance = 1)

for infection_start_point in [infection_node]:
     # Load network
     A = np.load('Diffusion/Results/Contagion/'+case+'_'+str(node_no)+'_network.npy')
     network = nx.from_numpy_matrix(A)
     # Load simulated dataset
     y_obs = np.load('Diffusion/Results/Contagion/'+case+'_yobs_'+str(infection_start_point)+'.npy')
     print('Infection Start point:' + str(infection_start_point))
     for ind in [0]:#range(len(y_obs)):
         kernel = Diffusion_kernel(network, mean = 0, cov = 1, seed_node = 1, seed = 1)
         #Redefine prior given observed dataset       
         prior = DiffusionPrior(np.array([1,1,0]), network, np.array(y_obs[ind][0][0]),[0.0,0.0,0.0],[1.0,1.0,1.0],seed = 1)
         #Initiate model with new prior
         model = ComplexContagion(prior, network, time_observed, theta = .3, gamma = .5, infection_start_point = [0], seed = 1)
         # Sampling using SABC       
         sampler_sabc = SABCDiffusion(model, dist_calc, kernel, backend, seed = 1)
         step, epsilon_sabc, n_samples, n_samples_per_param, ar_cutoff, beta, delta, v = 2, np.array([40.0]), 4, 1, 0.0001, 2, 0.2, 0.3
         start_time = time.time()
         journal_sabc = sampler_sabc.sample([y_obs[ind]], step, epsilon_sabc, n_samples, n_samples_per_param, beta, delta, v, ar_cutoff, resample=None, n_update=None, adaptcov=1, full_output=1)
         print("--- %s seconds ---" % (time.time() - start_time))
         journal_sabc.save('Diffusion/Results/Contagion/'+case+'_'+str(node_no)+'_joint_SABC_'+str(infection_start_point)+'_'+str(ind)+'.jrnl')
        
     del A, network, y_obs, kernel, prior, model, sampler_sabc, journal_sabc