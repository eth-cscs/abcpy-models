import numpy as np
from abcpy.distances import Distance
import networkx as nx

class SubsetDistance(Distance):
    def __init__(self, statistics, network, sp_distance=1):
        self.statistics_calc = statistics
        self.network = network
        self.sp_distance = sp_distance
        self.all_length = nx.all_pairs_shortest_path_length(network)
        self.max_length = max(self.all_length)

    def distance(self, d1, d2):
        S1n = self.statistics_calc.statistics(d1)
        S2n = self.statistics_calc.statistics(d2)

        result1, result2 = 0, 0
        for ind in range(len(S1n[0])):
            S1 = S1n[0][ind][0]
            S2 = S2n[0][ind][0]

            # Euclidean distance between proportion of infected
            result1 += pow(((len(S1) - len(S2)) / len(self.network.nodes())), 2)

            # Shortest path distance between subsets containing infected nodes
            if self.sp_distance == 1:
                result2 += self.local_distance(S1, S2, self.all_length)

        result = (np.sqrt(result1) + result2) / len(S1n[0])
        #print(result)
        return (result)

    def local_distance(self, S1, S2, all_length):

        # Take the disjoint parts of these two subsets
        S1_dj = np.setdiff1d(S1, S2)
        S2_dj = np.setdiff1d(S2, S1)

        result = 0
        for ind1 in range(len(S1_dj)):
            for ind2 in range(len(S2_dj)):
                result += all_length[S1_dj[ind1]][S2_dj[ind2]] / self.max_length
        if len(S1_dj) != 0 and len(S2_dj) != 0:
            result = result / (len(S1_dj) * len(S2_dj))
        elif len(S1_dj) * len(S2_dj) == 0 and len(S1_dj) + len(S2_dj) != 0:
            result = 1

        return (result)

    def dist_max(self):
        return self.max_length
