from abcpy.models import Model
import numpy as np
import networkx as nx
import copy
                       
class Diffusion_model_SI(Model):
    """Ecological model that describes the observed size of animal population over time 
    described in [1].
        
    [1] S. N. Wood. Statistical inference for noisy nonlinear ecological 
    dynamic systems. Nature, 466(7310):1102–1104, Aug. 2010.
    """

    def __init__(self, prior, network, time_observed, theta = None, infection_start_point = 0, seed = None):
        """
        Parameters
        ----------
        prior: abcpy.distributions.Distribution
            Prior distribution
        network: list 
            The fixed network structure on which deisease spreads
        time_observed: ndarray
            The timepoints where the observation were made
        theta: list or numpy.ndarray, optional       
            The parameter is a vector consisting of three numbers \
            :math:`\log r` (real number), :math:`\sigma` (positive real number, > 0), :math:`\phi` (positive real number > 0)
            If the parameter is ommitted, sampled from the prior.
        infection_start_point: integer
            The node where the infection starts (seed node)
        seed: int, optional
            Initial seed. The default value is generated randomly.             
        """    
        # test prior
        self.prior = prior
        self.node_no = len(network.nodes())
        self.network = network
        self.n_timestep = max(time_observed)
        self.time_observed = time_observed
        
        # set model parameters directly if specified
        if theta != None:
            if self.set_parameters(np.concatenate(([theta], [infection_start_point]))) == False:
                raise ValueError("The parameter values are out of the model parameter domain.")
        else:
            self.sample_from_prior()
        self.rng = np.random.RandomState(seed)
                    
    def sample_from_prior(self):
        sample = self.prior.sample(1).reshape(-1)
        if self.set_parameters(sample) == False:
            raise ValueError("Prior generates values that are out the model parameter domain.")                

    def simulate(self, n_simulate):
        diffusionstate_array = [None]*n_simulate
        # Initialize local parameters  
        for k in range(n_simulate):   
            # Initialize the time-series
            diffusionstate_array_tmp = []
            infected_nodes = list(self.infection_start_point)
            present_infected_nodes = copy.deepcopy(infected_nodes)
            if 0 in self.time_observed:                
                diffusionstate_array_tmp.append([present_infected_nodes])
            for ind_t in range(1,self.n_timestep):
                for ind_l in present_infected_nodes:
                    chosen_node_for_infection = self.rng.choice(self.network.neighbors(ind_l),1)[0]
                    if (chosen_node_for_infection in infected_nodes) == 0:
                        if self.rng.binomial(1,self.theta) == 1:
                            infected_nodes.append(chosen_node_for_infection)
                present_infected_nodes = copy.deepcopy(infected_nodes)
                if ind_t in self.time_observed:                                    
                    diffusionstate_array_tmp.append([present_infected_nodes])        
            diffusionstate_array[k] = diffusionstate_array_tmp
        # return an array of objects of list containing infected_nodes at the observed time-points
        return diffusionstate_array
    
    def get_parameters(self):
        return np.concatenate(([self.theta], self.infection_start_point))

    def set_parameters(self, theta):
        if isinstance(theta, (list,np.ndarray)):
                theta = np.array(theta)
        else:
            raise TypeError('The parameter value is not of allowed types')            
        #if theta.shape[0] > 1: return False
        if theta[0] < 0 or theta[0] > 1: return False
        self.theta = theta[0]
        if theta[1] < 0 or theta[1] > self.node_no-1: return False
        self.infection_start_point = theta[1:]
        return True
