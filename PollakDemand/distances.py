import math

import numpy as np
from abcpy.distances import Distance
from scipy.stats import entropy, wasserstein_distance


class KullbackLiebler(Distance):
    def __init__(self, statistics):
        self.statistics_cal = statistics
        
        self.s1 = None
        self.data_set = None
        self.dataSame = False
        
    def distance(self, d1, d2):
        """
        :param list d1: A list, containing a list describing the dataset
        :param list d2: A list, containing a list describing the data set
                
        :note: Hellinger operates on whole distribution not summary statistics.
        """
        if not isinstance(d1, list):
            raise TypeError('Data is not of allowed types')
        if not isinstance(d2, list):
            raise TypeError('Data is not of allowed types')
        
        #Check whether d1 is same as self.data_set
        if self.data_set is not None:
            if len(np.array(d1[0]).reshape(-1,)) ==1:
                self.data_set == d1
            else:
                self.dataSame = all([(np.array(self.data_set[i]) == np.array(d1[i])).all() for i in range(len(d1))])
        
        #Summary Statistics set to datasets
        if(self.s1 is None or self.dataSame is False):
            self.s1 = d1
            self.data_set = d1
        s2 = d2

        # Compute distnace between the statistics (whole dataset)
        dist = entropy(np.concatenate(self.s1,axis=0).reshape(-1,) , np.concatenate(s2,axis=0).reshape(-1,) )
        return dist
    
    def dist_max(self):
        return np.inf
    
class Hellinger(Distance):
    def __init__(self, statistics):
        self.statistics_cal = statistics
        
        self.s1 = None
        self.data_set = None
        self.dataSame = False
        
    def distance(self, d1, d2):
        """Parameters
            ----------
            d1, d2: list
                A list, containing a list describing the data set
            
        Note
            ----------
            Hellinger operates on whole distribution not summary statistics.
        """
        def hellinger_calculator(p, q):
            list_of_squares = []
            for p_i, q_i in zip(p, q):
                s = (math.sqrt(p_i) - math.sqrt(q_i)) ** 2
                list_of_squares.append(s)
            sosq = sum(list_of_squares)    
            return sosq / math.sqrt(2)

        if not isinstance(d1, list):
            raise TypeError('Data is not of allowed types')
        if not isinstance(d2, list):
            raise TypeError('Data is not of allowed types')
        
        #Check whether d1 is same as self.data_set
        if self.data_set is not None:
            if len(np.array(d1[0]).reshape(-1,)) ==1:
                self.data_set == d1
            else:
                self.dataSame = all([(np.array(self.data_set[i]) == np.array(d1[i])).all() for i in range(len(d1))])
        
        #Summary Statistics set to datasets
        if(self.s1 is None or self.dataSame is False):
            self.s1 = d1
            self.data_set = d1

        s2 = d2

        # Compute distance between the datasets
        dist = hellinger_calculator(np.concatenate(self.s1,axis=0).reshape(-1,).tolist(),np.concatenate(s2,axis=0).reshape(-1,).tolist()  ) 
        
        return dist
    
    def dist_max(self):
        return np.inf
    
class Wasserstein_distance(Distance):
    def __init__(self, statistics):
        self.statistics_cal = statistics
        
        self.s1 = None
        self.data_set = None
        self.dataSame = False
        
    def distance(self, d1, d2):
        """Parameters
            ----------
            d1, d2: list
                A list, containing a list describing the data set
            
            Note
              ----------
              Wasserstein operates on whole distribution not summary statistics.
        """
       
        if not isinstance(d1, list):
            raise TypeError('Data is not of allowed types')
        if not isinstance(d2, list):
            raise TypeError('Data is not of allowed types')
        
        #Check whether d1 is same as self.data_set
        if self.data_set is not None:
            if len(np.array(d1[0]).reshape(-1,)) ==1:
                self.data_set == d1
            else:
                self.dataSame = all([(np.array(self.data_set[i]) == np.array(d1[i])).all() for i in range(len(d1))])
        
        #Summary Statistics set to datasets
        if(self.s1 is None or self.dataSame is False):
            self.s1 = d1
            self.data_set = d1
        s2 = d2

        # Compute distance between the datasets
        dist = wasserstein_distance(np.concatenate(self.s1,axis=0).reshape(-1,),np.concatenate(s2,axis=0).reshape(-1,)  ) 
        
        return dist
    
    def dist_max(self):
        return np.inf
