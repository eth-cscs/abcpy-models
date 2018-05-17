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
    def __init__(self, network, sample, x_sn):
        self.network = network
        self.sample = sample
        self.unique_sn = np.unique(self.sample[:, 1])
        self.L = sample.shape[0]
        self.x_sn = x_sn

    def compute(self, x):
        result = 0
        for ind in self.unique_sn:
            marginal = sum(self.sample[:, 1] == ind) / self.L
            result += marginal * np.mean(
                np.sqrt(np.sum(pow(self.sample[self.sample[:, 1] == ind, 0] - x, 2))) \
                + nx.shortest_path_length(self.network, ind, self.x_sn))
        return (result)

    def minimize(self):
        Optim_sol = np.zeros(shape=(len(self.unique_sn),))
        Optim_val = np.zeros(shape=(len(self.unique_sn),))
        for ind in range(len(self.unique_sn)):
            self.x_sn = self.unique_sn[ind]
            res = minimize(self.compute, 0.5, tol=1e-6)
            Optim_sol[ind] = res.x
            Optim_val[ind] = self.compute(res.x)
        min_index = np.where(Optim_val == min(Optim_val))[0][0]
        Result = [Optim_sol[min_index], self.unique_sn[min_index]]
        return (Result)