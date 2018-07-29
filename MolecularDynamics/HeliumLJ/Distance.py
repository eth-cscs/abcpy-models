from abcpy.distances import Distance
from scipy import stats
import numpy as np

class KLdistance(Distance):
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
                kernel1 = stats.gaussian_kde(self.s1[ind1,:])
                kernel2 = stats.gaussian_kde(s2[ind2,:])
                xmin, xmax, gridsize = min(min(self.s1[ind1, :]), min(s2[ind2, :])), max(max(self.s1[ind1,:]),max(s2[ind2,:])), 100
                position = np.linspace(xmin, xmax, gridsize)
                dist[ind1, ind2] = stats.entropy(kernel1(position), kernel2(position))

        return dist.mean()

    def dist_max(self):
        return np.inf