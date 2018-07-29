import numpy as np
from abcpy.distances import Distance

class Abserror(Distance):
    """
    This class implements the Absolute error distance between two vectors.

    The maximum value of the distance is np.inf.
    """

    def __init__(self, statistics):
        self.statistics_calc = statistics

        # Since the observations do always stay the same, we can save the summary statistics of them and not recalculate it each time
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
        d1 = d1[0]
        d2 = d2[0]
        # Extract summary statistics from the dataset
        if (self.s1 is None):
            self.s1 = self.statistics_calc.statistics(d1)
            self.data_set = d1
        s2 = self.statistics_calc.statistics(d2)

        # compute distance between the statistics
        dist = np.zeros(shape=(self.s1.shape[0], s2.shape[0]))
        for ind1 in range(0, self.s1.shape[0]):
            for ind2 in range(0, s2.shape[0]):
                dist[ind1, ind2] = np.mean(np.abs(self.s1[ind1, :] - s2[ind2, :]))

        return dist.mean()

    def dist_max(self):
        return np.inf

