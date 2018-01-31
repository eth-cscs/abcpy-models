from abcpy.models import Model
import numpy as np
from abcpy.statistics import Statistics
from abcpy.distributions import Distribution
from scipy.stats import norm
from abcpy.distances import Distance
import networkx as nx
from abcpy.output import Journal
from scipy import optimize
import copy


class Network_generate:
    
    def __init__(self, node_no = 100, seed=None):
        """
        Defines the upper and lower bounds of a p-dimensional uniform Prior distribution in a closed interval.

        Parameters
        ----------
        node_no: integer 
            The number of nodes in the network
        seed: integer
            Initial seed for the random number generator
        """

        self.node_no = node_no
        self.rng = np.random.RandomState(seed)
        
    def linear(self):
        
        network = []
        
        for ind1 in np.arange(0,self.node_no-1):
            network.append([ind1,ind1+1])
        # Create a networkx Graph element
        G = nx.Graph()
        G.add_edges_from(network)            
        return(G)

        
    def lattice_2d(self,m):
         
        network = []
         
        for ind1 in range(m):
            for ind2 in range(int(self.node_no/m-1)):
                network.append([int((self.node_no/m)*ind1+ind2),int((self.node_no/m)*ind1+ind2+1)])
        
        for ind1 in range(m-1):
            for ind2 in range(int(self.node_no/m)):
                network.append([int((self.node_no/m)*ind1+ind2),int((self.node_no/m)*ind1+ind2+(self.node_no/m))])
        # Create a networkx Graph element
        G = nx.Graph()
        G.add_edges_from(network)                            
        return(G)
        
    def erdos_renyi(self, p = 0.3):
        
        # Generate a Erdo-renyi Graph
        G = nx.erdos_renyi_graph(self.node_no, p, seed=self.rng.randint(np.iinfo(np.int32).max)) 
        # Get the giant component
        giant = max(nx.connected_component_subgraphs(G), key=len)
        return(giant)
        
        
    def barabasi_albert(self, m = 2):
              
        return(nx.barabasi_albert_graph(self.node_no, m, seed=self.rng.randint(np.iinfo(np.int32).max)))

