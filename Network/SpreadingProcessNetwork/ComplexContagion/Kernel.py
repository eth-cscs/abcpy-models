import numpy as np
from abcpy.distributions import Distribution
import networkx as nx


class DiffusionKernel(Distribution):
    """
    This class implements a kernel for diffusion, which is 1-dimensional Normal distribution
    in first coordinate and second one is a discrete ditribution chosing connected nodes with 
    equal probability.
    """

    def __init__(self, network, mean=0, cov=1, seed_node=1, seed=None):
        """        
        Defines and mean and variance of normal distribution.

        Parameters
        ----------
        mean: numpy.ndarray
            mean
        var: numpy.ndarray
            variance
        network: list
            Structure of the network consisting of a list of connected edges
        seed: integer
            Initial seed for the random number generator

        """

        self.rng = np.random.RandomState(seed)
        self.mean = mean
        self.cov = cov
        self.seed_node = seed_node
        self.network = network
        self.all_degree = np.zeros(shape=(len(network.nodes()),))
        for ind in range(len(self.all_degree)):
            self.all_degree[ind] = network.degree(ind)

    def set_parameters(self, params):
        self.mean = params[0]
        self.cov = params[1]
        self.seed_node = params[2]

    def reseed(self, seed):
        self.rng.seed(seed)

    def sample(self, k):
        # samples = self.distribution.rvs(k).reshape(k,p)
        samples_1 = self.rng.multivariate_normal(self.mean, self.cov, k)
        nodes_proposed = self.network.neighbors(self.seed_node)
        weight = np.zeros(shape=(len(nodes_proposed),))
        for ind in range(len(nodes_proposed)):
            weight[ind] = 1 / self.all_degree[nodes_proposed[ind]]
        weight = weight / sum(weight)
        samples_2 = np.random.choice(nodes_proposed, k, p=weight)
        return np.transpose(np.reshape(np.append(samples_1, samples_2), (3, k)))

    def pdf(self, x):
        if x[2] in self.network.neighbors(self.seed_node):
            return multivariate_normal.pdf(x[0, 1], self.mean, self.cov) * (
            1 / len(self.network.neighbors(self.seed_node)))
        else:
            return 0
