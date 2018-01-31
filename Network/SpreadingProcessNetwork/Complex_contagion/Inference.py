import numpy as np
import networkx as nx
from Diffusion.Diffusion_contagion import Diffusion_model_complex_contagion, SABC_Diffusion, Diffusion_kernel, Diffusion_model_complex_contagion_theta
#from Diffusion_network_distance import Diffusion_Identity_statistics, SI_distance, SIR_distance
from Diffusion.Diffusion import Diffusion_prior, Network_generate, Diffusion_Identity_statistics, Subset_distance
import time
#import matplotlib.pyplot as plt
#from scipy.stats import gaussian_kde
from abcpy.distributions import Uniform, MultiNormal
from abcpy.inferences import SABC
#import matplotlib.pyplot as plt
#from scipy.stats import gaussian_kde

from abcpy.backends import BackendDummy
backend = BackendDummy()

#from abcpy.backend_mpi import BackendMPI as Backend
#backend = Backend()

# Different cases
case, node_no = 'ba', 100 # other cases 'ba', 'inrv', 'fb', 'tw', 'gog'
# Time observed
time_observed = np.arange(100, 200)
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
#theta, gamma, seed = 0.3, 0.5, 1
# Complex Contagion process simulation
#prior = Diffusion_prior(np.array([1,1,0]), network, network.nodes(),[0.0,0.0,0.0],[1.0,1.0,1.0],seed = 1) 
#for infection_start_point in range(node_no):   
#    model = Diffusion_model_complex_contagion(prior, network, time_observed, theta, gamma, [infection_start_point], seed)
#    y_obs = model.simulate(100)    
#    np.save('Diffusion/Results/Contagion/'+case+'_yobs_'+str(infection_start_point)+'.npy',y_obs)

# Inference
#==============================================================================
# 
# Import Statistics and Distance
stat_calc = Diffusion_Identity_statistics(degree=1, cross=0)
dist_calc = Subset_distance(stat_calc, network)


for infection_start_point in range(100):
     # Load network
     A = np.load('Diffusion/Results/Contagion/'+case+'_'+str(node_no)+'_network.npy')
     network = nx.from_numpy_matrix(A)
     # Load simulated dataset
     y_obs = np.load('Diffusion/Results/Contagion/'+case+'_yobs_'+str(infection_start_point)+'.npy')
     print('Infection Start point:' + str(infection_start_point))
     for ind in range(len(y_obs)):
         print(ind)        
         if len(y_obs[ind][0])==1:
             #Define kernle 
             kernel = MultiNormal(np.array([-13.0, .0, 7.0]), np.eye(3), seed=1)
             #Define prior
             prior = Uniform([0.0, 0.0],[1.0, 1.0],seed = 1)
             #Initiate Model
             model = Diffusion_model_complex_contagion_theta(prior, network, time_observed, theta = .3, gamma = .5, infection_start_point = y_obs[ind][0], seed = 1)
             # Sampling using SABC
             sampler_sabc = SABC(model, dist_calc, kernel, backend, seed = 1)
             step, epsilon_sabc, n_samples, n_samples_per_param, eps_percentile = 2, np.array([40.0]), 10, 1, 10
             start_time = time.time()
             journal_sabc = sampler_sabc.sample([y_obs[ind]], step, epsilon_sabc, n_samples, n_samples_per_param, eps_percentile)
             print("--- %s seconds ---" % (time.time() - start_time))
             journal_sabc.save('Diffusion/Results/Contagion/'+case+'_'+str(node_no)+'_joint_SABC_'+str(infection_start_point)+'_'+str(ind)+'.jrnl')
         else:    
             kernel = Diffusion_kernel(network, mean = 0, cov = 1, seed_node = 1, seed = 1)
             #Redefine prior given observed dataset       
             prior = Diffusion_prior(np.array([1,1,0]), network, np.array(y_obs[0][0]),[0.0,0.0,0.0],[1.0,1.0,1.0],seed = 1)
             #Initiate model with new prior
             model = Diffusion_model_complex_contagion(prior, network, time_observed, theta = .3, gamma = .5, infection_start_point = [0], seed = 1)    
             # Sampling using SABC       
             sampler_sabc = SABC_Diffusion(model, dist_calc, kernel, backend, seed = 1)
             step, epsilon_sabc, n_samples, n_samples_per_param, eps_percentile = 2, np.array([40.0]), 10, 1, 10
             start_time = time.time()
             journal_sabc = sampler_sabc.sample([y_obs[ind]], step, epsilon_sabc, n_samples, n_samples_per_param, eps_percentile)
             print("--- %s seconds ---" % (time.time() - start_time))
             journal_sabc.save('Diffusion/Results/Contagion/'+case+'_'+str(node_no)+'_joint_SABC_'+str(infection_start_point)+'_'+str(ind)+'.jrnl')
        
     del A, network, y_obs, kernel, prior, model, sampler_sabc, journal_sabc
 
#==============================================================================
