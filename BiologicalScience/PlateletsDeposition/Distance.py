class DepositionDistance(Distance):
    """
    This class implements the Euclidean distance between two vectors.

    The maximum value of the distance is np.inf.
    """
    
    def __init__(self, statistics):
        self.statistics_calc = statistics

        
    def distance(self, d1, d2):
        if len(d1) != len(d2):
            raise BaseException("Input data sets have different sizes: {} vs {}".format(len(d1), len(d2)))

        s1 = self.statistics_calc.statistics(d1)
        s2 = self.statistics_calc.statistics(d2)
        
        result = 0
        for ind in range(4):
            result += 0.25*(1 - np.exp(-self._bhattacharya_coefficient(s1[0,ind], s1[0,ind+4], s2[0,ind], s2[0,ind+4])))
        result = 0.5*(result) + 0.5*np.sqrt(np.mean(pow((s1[0,8:]-s2[0,8:]),2)))
               
        return result

    
    def dist_max(self):
        return np.inf

    def _bhattacharya_coefficient(self, mean_x, var_x, mean_y, var_y):
        
        return 0.25*np.log(0.25*(var_x/var_y + var_y/var_x + 2)) + 0.25*(pow((mean_x-mean_y),2)/(var_x + var_y)) 

class DepositionRelativeErrorDistance(Distance):
    """
    This class implements the Euclidean distance between two vectors.

    The maximum value of the distance is np.inf.
    """
    
    def __init__(self, statistics):
        self.statistics_calc = statistics

        
    def distance(self, d1, d2):
        if len(d1) != len(d2):
            raise BaseException("Input data sets have different sizes: {} vs {}".format(len(d1), len(d2)))

        s1 = self.statistics_calc.statistics(d1)
        s2 = self.statistics_calc.statistics(d2)
        
        error = np.zeros(shape=(s1.shape[1],1))

        for ind1 in range(s1.shape[1]):
            for ind2 in range(s1.shape[0]):
                if s1[ind2,ind1] != 0:
                    error[ind1] += np.abs(s1[ind2,ind1]-s2[ind2,ind1])/np.abs(max(s1[ind2,ind1],s2[ind2,ind1]))
                else:
                    error[ind1] += np.abs(s1[ind2,ind1]-s2[ind2,ind1])
        
        error *= 1/5
        return np.mean(error)
    
    def dist_max(self):
        return 1



class DepositionDistanceCombined(Distance):
    """
    This class implements the Euclidean distance between two vectors.

    The maximum value of the distance is np.inf.
    """
    
    def __init__(self, statistics):
        self.statistics_calc = statistics

        
    def distance(self, d1, d2):
        #if len(d1) != len(d2):
        #    raise BaseException("Input data sets have different sizes: {} vs {}".format(len(d1), len(d2)))

        s1_all = self.statistics_calc.statistics(d1)
        s2_all = self.statistics_calc.statistics(d2)
        
        s1 = s1_all[0]
        s2 = s2_all[0]        
        
        result = 0
        for ind in range(4):
            result += 0.25*(1 - np.exp(-self._bhattacharya_coefficient(s1[0,ind], s1[0,ind+4], s2[0,ind], s2[0,ind+4])))
        result = 0.5*(result) + 0.5*np.sqrt(np.mean(pow((s1[0,8:]-s2[0,8:]),2)))
        
        s1 = s1_all[1]
        s2 = s2_all[1]
        
        error = np.zeros(shape=(s1.shape[1],1))

        for ind1 in range(s1.shape[1]):
            for ind2 in range(s1.shape[0]):
                if s1[ind2,ind1] != 0:
                    error[ind1] += np.abs(s1[ind2,ind1]-s2[ind2,ind1])/np.abs(max(s1[ind2,ind1],s2[ind2,ind1]))
                else:
                    error[ind1] += np.abs(s1[ind2,ind1]-s2[ind2,ind1])
        
        error *= 1/5
        
        return np.mean(error) + result

    
    def dist_max(self):
        return np.inf

    def _bhattacharya_coefficient(self, mean_x, var_x, mean_y, var_y):
        
        return 0.25*np.log(0.25*(var_x/var_y + var_y/var_x + 2)) + 0.25*(pow((mean_x-mean_y),2)/(var_x + var_y)) 
