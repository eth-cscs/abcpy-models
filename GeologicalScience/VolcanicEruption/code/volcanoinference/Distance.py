import numpy as np
from fileio import *
from abcpy.distances import Distance

# different types of distance computation between extracted vectors
def meansquarederror(measures, compute):
    return mean((measures-compute)**2)

def meanerror(measures, compute):
    return mean((measures-compute))

def meanr(measures, compute):
    ss_res = dot((measures - compute),(measures - compute))
    ymean = mean(measures)
    ss_tot = dot((measures-ymean),(measures-ymean))
    return 1.0-ss_res/ss_tot

# http://www.netlib.org/lapack/lug/node75.html
def relativeerror(measures, compute):
    return linalg.norm(compute-measures)/linalg.norm(measures)

# https://en.wikipedia.org/wiki/Mean_absolute_percentage_error
def meanrelativeerror(measures, compute):
    return mean(abs(array(measures-compute))/array(measures))

class DistanceVolcano(Distance):
    """
    This class implements the Euclidean distance between two vectors.

    The maximum value of the distance is np.inf.
    """

    def __init__(self, statistics):
        self.statistics_calc = statistics

        # Since the observations do always stay the same, we can save the
        #  summary statistics of them and not recalculate it each time
        self.s1 = None
        self.data_set = None

    def distance(self, d1, d2):
        """Calculates the distance between two datasets.

        Parameters
        ----------
        d1, d2: list
            A list, containing a list describing the data set
        """
        if not isinstance(d1, list):
            raise TypeError('Data is not of allowed types')
        if not isinstance(d2, list):
            raise TypeError('Data is not of allowed types')

        # Extract summary statistics from the dataset
        if (self.s1 is None or self.data_set != d1):
            self.s1 = self.statistics_calc.statistics(d1)
            self.data_set = d1
        s2 = self.statistics_calc.statistics(d2)

        # compute distance between the statistics
        dist = np.zeros(shape=(self.s1.shape[0], s2.shape[0]))
        for ind1 in range(0, self.s1.shape[0]):
            for ind2 in range(0, s2.shape[0]):
                dist[ind1, ind2] = meansquarederror(self.s1[ind1, :], s2[ind2, :])

        return dist.mean()

    def dist_max(self):
        return np.inf

