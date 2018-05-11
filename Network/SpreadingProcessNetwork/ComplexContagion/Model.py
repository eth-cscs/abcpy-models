import copy
import numpy as np
import networkx as nx
from scipy import optimize
from abcpy.models import Model

class ComplexContagion(Model):
    """Complex Contagion model for fake news spread on facebook social network.
    """

    def __init__(self, prior, network, time_observed, theta=None, gamma=None, infection_start_point=[0], seed=None):
        """
        Parameters
        ----------
        prior: abcpy.distributions.Distribution
            Prior distribution
        network: list
            The fixed network structure on which deisease spreads
        time_observed: ndarray
            The timepoints where the observation were made
        theta: float
            Rate of exposure
        gamma: float
            Relative threshold
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
            if self.set_parameters(np.concatenate(([theta, gamma], infection_start_point))) == False:
                raise ValueError("The parameter values are out of the model parameter domain.")
        else:
            self.sample_from_prior()
        self.rng = np.random.RandomState(seed)

    def sample_from_prior(self):
        sample = self.prior.sample(1).reshape(-1)
        if self.set_parameters(sample) == False:
            raise ValueError("Prior generates values that are out the model parameter domain.")

    def simulate(self, n_simulate):
        diffusionstate_array = [None] * n_simulate
        for k in range(n_simulate):
            # All neighbors
            neighbors = self.network.neighbors(self.infection_start_point[0])
            # We choose gamma proportion of neighbors as native adaptors
            chosen_neighbors = list(self.rng.choice(neighbors, int(len(neighbors) * self.gamma)))
            chosen_neighbors.append(int(self.infection_start_point[0]))
            # Infected nodes to start with
            infected_nodes = list(chosen_neighbors)
            # infected_nodes = self.infection_start_point
            exposed_nodes, count_of_exposure, diffusionstate_array_tmp = [], [], []
            present_infected_nodes = copy.deepcopy(infected_nodes)
            if 0 in self.time_observed:
                diffusionstate_array_tmp.append(
                    [present_infected_nodes, copy.deepcopy(exposed_nodes), copy.deepcopy(count_of_exposure)])
            # Loop over the time-steps
            for ind_t in range(1, self.n_timestep):
                # Loop over all the infected_nodes
                for ind_infected in present_infected_nodes:
                    chosen_node_for_exposure = self.rng.choice(self.network.neighbors(ind_infected), 1)[0]
                    if (chosen_node_for_exposure in infected_nodes) == 0:
                        # The chosen node got exposed
                        if self.rng.binomial(1, self.theta) == 1:
                            if (chosen_node_for_exposure in exposed_nodes):
                                count_of_exposure[exposed_nodes.index(chosen_node_for_exposure)] \
                                    [len(set(infected_nodes).intersection(
                                    self.network.neighbors(chosen_node_for_exposure))) - 1] += 1
                            else:
                                exposed_nodes.append(chosen_node_for_exposure)
                                count_of_exposure.append(
                                    list(np.zeros(shape=(len(self.network.neighbors(chosen_node_for_exposure)),))))
                                count_of_exposure[exposed_nodes.index(chosen_node_for_exposure)] \
                                    [len(set(infected_nodes).intersection(
                                    self.network.neighbors(chosen_node_for_exposure))) - 1] += 1
                                ## Check whether the exposed node got infected
                            # Count of exposure for the exposed node
                            L = count_of_exposure[exposed_nodes.index(chosen_node_for_exposure)]
                            if self.decision_adoption(L, len(L)):
                                infected_nodes.append(chosen_node_for_exposure)
                                # del count_of_exposure[exposed_nodes.index(chosen_node_for_exposure)]
                                # del exposed_nodes[exposed_nodes.index(chosen_node_for_exposure)]
                present_infected_nodes = copy.deepcopy(infected_nodes)
                if ind_t in self.time_observed:
                    diffusionstate_array_tmp.append(
                        [present_infected_nodes, copy.deepcopy(exposed_nodes), copy.deepcopy(count_of_exposure)])
            diffusionstate_array[k] = diffusionstate_array_tmp

        return diffusionstate_array

    def decision_adoption(self, LL, F):
        probability = 1
        L = LL[:max(np.where(np.array(LL) > 0)[0]) + 1]
        for ind in range(len(L) - 1):
            probability *= pow((1 - self.likelihood_adoption(ind, F)), L[ind])
        probability *= pow((1 - self.likelihood_adoption(len(L), F)), L[-1] - 1) * self.likelihood_adoption(len(L), F)

        return self.rng.binomial(1, probability)

    def likelihood_adoption(self, k, F):
        epsilon_high, epsilon_low, g = 0.25, 0.001, 1
        prob = epsilon_low + ((epsilon_high - epsilon_low) / (1 + np.exp(-g * F * (k / F - self.gamma))))
        return prob

    def get_parameters(self):
        return np.concatenate(([self.theta, self.gamma], self.infection_start_point))

    def set_parameters(self, theta):
        if isinstance(theta, (list, np.ndarray)):
            theta = np.array(theta)
        else:
            raise TypeError('The parameter value is not of allowed types')
            # if theta.shape[0] > 1: return False
        if theta[0] < 0 or theta[0] > 1: return False
        self.theta = theta[0]
        if theta[1] < 0 or theta[1] > 1: return False
        self.gamma = theta[1]
        if theta[2] < 0 or theta[2] > self.node_no - 1: return False
        self.infection_start_point = [int(i) for i in list(theta[2:])]
        return True