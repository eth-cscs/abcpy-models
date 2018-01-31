import numpy as np
from abcpy.distances import Distance
import networkx as nx

class Subset_distance(Distance):
    
    def __init__(self, statistics, network, sp_distance = 1):
        self.statistics_calc = statistics
        self.network = network
        self.sp_distance = sp_distance
        self.all_length = nx.all_pairs_shortest_path_length(network)
        self.max_length = max(self.all_length)
        
    def distance(self, d1, d2):
        S1n = self.statistics_calc.statistics(d1)
        S2n = self.statistics_calc.statistics(d2)

        result1, result2, result3, result4, result5 = 0, 0, 0, 0, 0
        for ind in range(len(S1n[0])):

            # Euclidean distance between proportion of infected
            result1 += pow(((len(S1n[0][ind][0])-len(S2n[0][ind][0]))/len(self.network.nodes())),2)

            # Euclidean distance between proportion of exposed
            result2 += pow(((len(S1n[0][ind][1])-len(S2n[0][ind][1]))/len(self.network.nodes())),2)
            
            # Shortest path distance between subsets containing infected nodes
            # Shortest path distance between subsets containing exposed nodes
            if self.sp_distance == 1:
                result3 += self.local_distance(S1n[0][ind][0], S2n[0][ind][0], self.all_length)
                result4 += self.local_distance(list(set(S1n[0][ind][1])-set(S1n[0][ind][0])), list(set(S2n[0][ind][1])-set(S2n[0][ind][0])), self.all_length)
                
        result1 = np.sqrt(result1) 
        result2 = np.sqrt(result2)
        
        # Euclidean distance between the change in proportion of exposure
        result5 = self.local_exposure_rate_estimate(S1n[0], S2n[0])   
        
        result = (result1 + result2 + result3 + result4 + result5)/len(S1n[0])
        
        #print(result)
        return(result)

    def local_distance(self, S1, S2, all_length):
        
        # Take the disjoint parts of these two subsets
        S1_dj = np.setdiff1d(S1,S2)
        S2_dj = np.setdiff1d(S2,S1)
        
        result = 0
        for ind1 in range(len(S1_dj)):
            for ind2 in range(len(S2_dj)):
                result +=  all_length[S1_dj[ind1]][S2_dj[ind2]]/self.max_length
        if len(S1_dj) != 0 and len(S2_dj) != 0:
            result = result/(len(S1_dj)*len(S2_dj))
        elif len(S1_dj)*len(S2_dj) == 0:
            result = 1
        
        return(result)
        
    def local_exposure_rate_estimate(self, S1, S2):
        
        a = np.zeros(shape=(len(S1),1))
        b = np.zeros(shape=(len(S2),1))
        for ind in range(len(S1)):
            for ind1 in range(len(S1[ind][2])):
                a[ind] += sum(S1[ind][2][ind1])
            for ind1 in range(len(S2[ind][2])):
                b[ind] += sum(S2[ind][2][ind1])  
        result = 0
        for ind in np.arange(1,len(S1)):
            #result += pow(((a[ind]-a[ind-1])/len(S1[ind][0])) - ((b[ind]-b[ind-1])/len(S2[ind][0])),2)
            result += pow(((a[ind]-a[ind-1])-(b[ind]-b[ind-1]))/len(self.network.nodes()),2)
        return(np.sqrt(result))
        
    def dist_max(self):
        return self.max_length              
