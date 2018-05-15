import networkx as nx
import numpy as np
from scipy.optimize import minimize

class BayesEstimateComplexContagion:
    def __init__(self, network, sample, x_sn):
        self.network = network
        self.sample = sample
        self.unique_sn = np.unique(self.sample[:, 2])
        self.L = sample.shape[0]
        self.x_sn = x_sn

    def compute(self, x):
        result = 0
        for ind in self.unique_sn:
            marginal = sum(self.sample[:, 2] == ind) / self.L
            result += marginal * np.mean(
                np.sqrt(np.sum(pow(self.sample[self.sample[:, 2] == ind, :2] - x[:2], 2), axis=1)) \
                + nx.shortest_path_length(self.network, ind, self.x_sn))
        return (result)

    def minimize(self):
        Optim_sol = np.zeros(shape=(len(self.unique_sn), 2))
        Optim_val = np.zeros(shape=(len(self.unique_sn), 1))
        for ind in range(len(self.unique_sn)):
            self.x_sn = self.unique_sn[ind]
            res = minimize(self.compute, [0.5, 0.5], bounds=((0, 1), (0, 1)), tol=1e-6)
            Optim_sol[ind, :] = res.x
            Optim_val[ind] = self.compute(res.x)
        Result = np.concatenate((Optim_sol[np.where(Optim_val == min(Optim_val))[0][0], :], \
                                 self.unique_sn[np.where(Optim_val == min(Optim_val))[0][0]].reshape(1, )))
        return (Result)

class BayesEstimateSimpleContagion:
    def __init__(self, network, sample):
        self.network = network
        self.sample = sample

    def compute(self):
        a = np.unique(self.sample[:,1])
        post_prob = np.zeros(shape=(len(a),1))
        for ind1 in range(len(a)):
            post_prob[ind1] = np.sum(self.sample[:,1]==a[ind1])
        post_prob = post_prob/sum(post_prob)
        post_dist = np.zeros(shape=(len(a),1))
        for ind1 in range(len(a)):
            for ind2 in range(len(a)):
                post_dist[ind1] += post_prob[ind2]*nx.shortest_path_length(self.network, a[ind2], a[ind1])
        result = np.array([np.mean(self.sample[:,0]), a[np.where(post_dist == min(post_dist))[0][0]]])
        return (result)
