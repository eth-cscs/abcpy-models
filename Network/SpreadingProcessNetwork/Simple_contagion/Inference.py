import networkx as nx
import numpy as np 
from Diffusion.Diffusion_SC import Network_generate, Diffusion_prior, Diffusion_model_SI, Diffusion_kernel, SABC_Diffusion
import time

from abcpy.backends import BackendDummy
backend = BackendDummy()


# Different cases
case, infection_node = 'fb', 2000 # other cases 'ba', 'inrv', 'fb', 'tw', 'gog'
node_no, theta = 4039, 0.3
# Time observed
time_observed = np.arange(20, 70+1)
#time_observed = np.array([15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,27, 28, 29, 30])
# Load network
A = np.load('Diffusion/Results/SC/'+case+'_'+str(node_no)+'_network.npy')
network = nx.from_numpy_matrix(A)

# Generate Network
#==============================================================================
#generator = Network_generate(node_no,1)
#network = generator.barabasi_albert(4)
#network.nodes()
#A = nx.to_numpy_matrix(network)
#np.save('Results/All/'+case+'_'+str(node_no)+'_network.npy',A)
# 
#==============================================================================
#==============================================================================
# Indian rural village rading and saving the network
#Real example
#path = "Diffusion/Data/datav4.0/Data/1. Network Data/Adjacency Matrices/"
#Gs = []
#for k in range(1, 11):
#    A = np.genfromtxt(path + "adj_allVillageRelationships_vilno_" + str(k) + ".csv", delimiter=',')
#    G = nx.to_networkx_graph(A)
#    Gs.append(G)
#network = G
#A = nx.to_numpy_matrix(network)
#node_no = len(network.nodes())
#np.save('Diffusion/Results/SC/inrv_'+str(node_no)+'_network.npy',A)
#==============================================================================


# Generate dataset
#==============================================================================
theta, seed = 0.3, 1
prior = Diffusion_prior(np.array([1,0]), network, network.nodes(),[0.0,0.0],[1.0,1.0],seed = 1)
for infection_start_point in [infection_node]:   
    model = Diffusion_model_SI(prior, network, time_observed, theta, infection_start_point, seed)
    y_obs = model.simulate(100)    
    np.save('Diffusion/Results/SC/'+case+'_'+str(node_no)+'_yobs_'+str(infection_start_point)+'.npy',y_obs)
# 
#==============================================================================

# Inference
#==============================================================================
# 
# Import Statistics and Distance
from Diffusion.Diffusion_SC import Diffusion_Identity_statistics
from Diffusion.Diffusion_SC import Subset_distance
stat_calc = Diffusion_Identity_statistics(degree=1, cross=0)
dist_calc = Subset_distance(stat_calc, network)
# 
# Define model to simulate observed dataset
theta, seed = 0.3, 1
for infection_start_point in [infection_node]:
    # Load network
    A = np.load('Diffusion/Results/SC/'+case+'_'+str(node_no)+'_network.npy')
    network = nx.from_numpy_matrix(A)
    # Load simulated dataset
    y_obs = np.load('Diffusion/Results/SC/'+case+'_'+str(node_no)+'_yobs_'+str(infection_start_point)+'.npy')
    print('Infection Start point:' + str(infection_start_point))
    for ind in range(len(y_obs)):
        print(ind)           
        kernel = Diffusion_kernel(network, mean = 0, var = 1, seed_node = 1, seed = 1)
        #Redefine prior given observed dataset       
        prior = Diffusion_prior(np.array([1,0]), network, np.array(y_obs[ind][0][0]),[0.0,0.0],[1.0,1.0],seed = 1)
        #Initiate model with new prior
        model = Diffusion_model_SI(prior, network, time_observed, theta, infection_start_point, seed = None)    
        # Sampling using SABC       
        sampler_sabc = SABC_Diffusion(model, dist_calc, kernel, backend, seed = 1)
        step, epsilon_sabc, n_samples, n_samples_per_param, ar_cutoff, beta, delta, v = 2, np.array([40.0]), 4, 1, 0.0001, 2, 0.2, 0.3
        start_time = time.time()
        journal_sabc = sampler_sabc.sample([y_obs[ind]], step, epsilon_sabc, n_samples, n_samples_per_param, beta, delta, v, ar_cutoff, resample=None, n_update=None, adaptcov=1, full_output=1)
        print("--- %s seconds ---" % (time.time() - start_time))
        journal_sabc.save('Diffusion/Results/SC/'+case+'_'+str(node_no)+'_joint_SABC_'+str(infection_start_point)+'_'+str(ind)+'.jrnl')
    
    del A, network, y_obs, kernel, prior, model, sampler_sabc, journal_sabc
