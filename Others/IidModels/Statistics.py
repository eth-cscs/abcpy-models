import numpy as np
from abcpy.statistics import Statistics


class BetaStatistics(Statistics):
    """
    """

    def __init__(self, previous_statistics=None):

        self.previous_statistics = previous_statistics

    def statistics(self, data):

        # pipeline: first call the previous statistics:
        if self.previous_statistics is not None:
            data = self.previous_statistics.statistics(data)
        # the first of the statistics need to take list as input, in order to match the API. Then actually the
        # transformations work on np.arrays. In fact the first statistic transforms the list to array. Therefore, the
        # following code needs to be called only if the self statistic is the first, i.e. it does not have a
        # previous_statistic element.
        else:
            data = self._check_and_transform_input(data)

        return np.concatenate(
            (np.mean(np.log(data), axis=1).reshape(-1, 1), np.mean(np.log(1 - data), axis=1).reshape(-1, 1)), axis=1)


class GammaStatistics(Statistics):
    """
    """

    def __init__(self, previous_statistics=None):

        self.previous_statistics = previous_statistics

    def statistics(self, data):

        # pipeline: first call the previous statistics:
        if self.previous_statistics is not None:
            data = self.previous_statistics.statistics(data)
        # the first of the statistics need to take list as input, in order to match the API. Then actually the
        # transformations work on np.arrays. In fact the first statistic transforms the list to array. Therefore, the
        # following code needs to be called only if the self statistic is the first, i.e. it does not have a
        # previous_statistic element.
        else:
            data = self._check_and_transform_input(data)

        return np.concatenate((np.mean(np.log(data), axis=1).reshape(-1, 1), np.mean(data, axis=1).reshape(-1, 1)),
                              axis=1)


class GaussianStatistics(Statistics):
    """
    """

    def __init__(self, previous_statistics=None):

        self.previous_statistics = previous_statistics

    def statistics(self, data):

        # pipeline: first call the previous statistics:
        if self.previous_statistics is not None:
            data = self.previous_statistics.statistics(data)
        # the first of the statistics need to take list as input, in order to match the API. Then actually the
        # transformations work on np.arrays. In fact the first statistic transforms the list to array. Therefore, the
        # following code needs to be called only if the self statistic is the first, i.e. it does not have a
        # previous_statistic element.
        else:
            data = self._check_and_transform_input(data)

        return np.concatenate((np.mean(data, axis=1).reshape(-1, 1),
                               np.mean(data ** 2, axis=1).reshape(-1, 1)), axis=1)
