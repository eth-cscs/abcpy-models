import warnings

import numpy as np
from abcpy.statistics import Statistics


class Autocorrelation(Statistics):
    # another possible implementation is using the function at the following link, however it is slower:
    # https://www.statsmodels.org/devel/generated/statsmodels.tsa.stattools.acovf.html

    def __init__(self, order):
        self.order = order

    def statistics(self, data):

        if isinstance(data, list):
            if np.array(data).shape == (len(data),):
                if len(data) == 1:
                    data = np.array(data).reshape(1, 1)
                data = np.array(data).reshape(len(data), 1)
            else:
                data = np.concatenate(data).reshape(len(data), -1)
        else:
            raise TypeError('Input data should be of type list, but found type {}'.format(type(data)))

        results = np.zeros((data.shape[0], self.order))

        data_size = np.int(data.shape[1] / 2)

        for sample in range(data.shape[0]):
            results[sample] = np.correlate(data[sample], data[sample], mode="same")[
                              data_size + 1: data_size + self.order + 1]

        if np.any(np.isinf(results)):
            warnings.warn("Infinity in computation of Autocorrelation statistics", RuntimeWarning)
        if np.any(np.isnan(results)):
            warnings.warn("nan in computation of Autocorrelation statistics", RuntimeWarning)

        results = results / (data.shape[1] - np.arange(1, self.order + 1))
        # rescale the different correlations according to their range
        if np.any(np.isinf(results)):
            warnings.warn("Infinity in computation of Autocorrelation statistics", RuntimeWarning)
        if np.any(np.isnan(results)):
            warnings.warn("nan in computation of Autocorrelation statistics", RuntimeWarning)
        return results
