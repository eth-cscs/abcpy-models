# from abcpy.distances import Distance
# import numpy as np
# import string
# import ot
#
from abcpy.distances import Distance
import numpy as np

class WeightedEuclidean(Distance):
    """
    This class implements the Euclidean distance between two vectors.

    The maximum value of the distance is np.inf.
    """
    def __init__(self, statistics, weight):
        """
        Parameters
        ----------
        statistics_calc : abcpy.statistics.Statistics
            Statistics extractor object that conforms to the Statistics class.
        """
        self.weight = weight
        super(WeightedEuclidean, self).__init__(statistics)

    def distance(self, d1, d2):
        """Calculates the distance between two datasets, by computing Euclidean distance between each element of d1 and
        d2 and taking their average.

        Parameters
        ----------
        d1: Python list
            Contains n1 data points.
        d2: Python list
            Contains n2 data points.

        Returns
        -------
        numpy.float
            The distance between the two input data sets.
        """
        s1, s2 = self._calculate_summary_stat(d1, d2)

        # compute distance between the statistics
        dist = np.zeros(shape=(s1.shape[0], s2.shape[0]))
        for ind1 in range(0, s1.shape[0]):
            for ind2 in range(0, s2.shape[0]):
                dist[ind1, ind2] = self._list_diff(s1[ind1, :], s2[ind2, :])
        return dist.mean()

    def dist_max(self):
        """
        Returns
        -------
        numpy.float
            The maximal possible value of the desired distance function.
        """
        return 84

    def _list_diff(self, x, y):
        print(x, y)
        return np.sqrt(sum(self.weight * pow(x - y, 2)))