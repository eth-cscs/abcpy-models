import warnings

import numpy as np
import scipy.stats
from abcpy.statistics import Statistics


class LorenzLargerStatistics(Statistics):

    def __init__(self, moments=3, auto_covariances=5, covariances=5, cross_covariances=5, degree=1, cross=False):
        """We compute here moments, covariances, cross covariances, autocovariances."""
        self.number_moments = moments
        self.number_auto_covariances = auto_covariances
        self.number_covariances = covariances
        self.number_cross_covariances = cross_covariances
        self.degree = degree
        self.cross = cross

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
        # Extract Summary Statistics; this here assumes explicitly we are using 40 variables.
        num_element, timestep = len(data), int(data[0].shape[0] / 40)
        result = np.zeros(shape=(num_element, self.number_moments + self.number_auto_covariances +
                                 self.number_covariances + 2 * self.number_cross_covariances))
        # Compute statistics
        if np.any(np.isinf(data)):
            warnings.warn("Infinity in the output of Lorenz95 model", RuntimeWarning)
        if np.any(np.isnan(data)):
            warnings.warn("nan in the output of Lorenz95 model", RuntimeWarning)

        for ind_element in range(0, num_element):
            # First convert the vector to the 40 dimensional timeseries
            data_ind_element = data[ind_element].reshape(40, timestep)
            # compute 10 central moments:
            moments = np.zeros(self.number_moments)
            moments[0] = np.mean(np.mean(data_ind_element, 1))
            for i in range(1, self.number_moments):
                moments[i] = np.mean(scipy.stats.moment(data_ind_element, axis=1, moment=i + 1))
                # second centered moment corresponds to variance

            # auto-covariances: covariance between one timeseries and itself with some lag
            # Here we repeat computation, as we could return the auto-covariances for all lags with
            # one single function call. It does not matter too much for now.
            auto_covariances = np.zeros(self.number_auto_covariances)
            for ind in range(0, data_ind_element.shape[0]):  # loop over the different components of timeseries.
                for i in range(self.number_auto_covariances):
                    auto_covariances[i] += self._covariance_lag(data_ind_element[ind, :],
                                                                data_ind_element[ind, :], lag=i + 1)
            auto_covariances /= data_ind_element.shape[0]

            # covariances between different timeseries with lag 0 and spatial distance up to 5, with cyclic
            # boundary conditions.
            covariances_neighbors = np.zeros(self.number_covariances)
            for ind in range(0, data_ind_element.shape[0]):  # loop over the different components of timeseries.
                for i in range(self.number_covariances):
                    # the % is for the cyclic conditions
                    covariances_neighbors[i] += self._covariance_lag(
                        data_ind_element[ind, :], data_ind_element[(ind + i + 1) % data_ind_element.shape[0], :], lag=0)
            covariances_neighbors /= data_ind_element.shape[0]

            # we now extract cross-covariances with both time lag and spatial distance: consider timeseries at distance
            # i and use lag i, by keeping them to be the same:
            cross_covariances_neighbors_right = np.zeros(self.number_cross_covariances)
            cross_covariances_neighbors_left = np.zeros(self.number_cross_covariances)
            for ind in range(0, data_ind_element.shape[0]):  # loop over the different components of timeseries.
                for i in range(self.number_cross_covariances):
                    # the % is for the cyclic conditions
                    cross_covariances_neighbors_right[i] += self._covariance_lag(
                        data_ind_element[ind, :], data_ind_element[(ind + i + 1) % data_ind_element.shape[0], :],
                        lag=i + 1)
                    cross_covariances_neighbors_left[i] += self._covariance_lag(
                        data_ind_element[(ind + i + 1) % data_ind_element.shape[0], :], data_ind_element[ind, :],
                        lag=i + 1)
            cross_covariances_neighbors_right /= data_ind_element.shape[0]
            cross_covariances_neighbors_left /= data_ind_element.shape[0]

            result[ind_element] = np.concatenate((moments, auto_covariances, covariances_neighbors,
                                                  cross_covariances_neighbors_right, cross_covariances_neighbors_left))

            for i, stat in enumerate(result[ind_element]):
                if np.any(np.isinf(stat)):
                    warnings.warn("Infinity in stat number {}".format(i + 1), RuntimeWarning)
                if np.any(np.isnan(stat)):
                    warnings.warn("nan in stat number {}".format(i + 1), RuntimeWarning)
            # Expand the data with polynomial expansion
        result = self._polynomial_expansion(result)
        return np.array(result)

    def _covariance_lag(self, x, y, lag=0):
        """ Computes cross-covariance between x and y. Using np.correlate is faster than using the for loop over
        elements. This function can be used to compute all kinds of covariances/autocovariancess/cross-covariances
        needed.

        Parameters
        ----------
        x: numpy.ndarray
            Vector of real numbers.
        y: numpy.ndarray
            Vector of real numbers.

        Returns
        -------
        numpy.ndarray
            Cross-covariance calculated between x and y.
        """
        assert x.shape[0] == y.shape[0]  # we consider case in which the timeseries all have same length
        # data_size = np.int(x.shape[0] / 2)  # it may work only for timeseries with even number of timesteps.
        # the proper definition of covariance also subtracts the means of the two series that you are considering.
        if lag == 0:
            return np.correlate(x, y, mode="valid")[0] / x.shape[0] - np.mean(x) * np.mean(y)
        else:
            return np.correlate(x[lag:], y[:-lag], mode="valid")[0] / (x.shape[0] - lag) - np.mean(x[lag:]) * np.mean(
                y[:-lag])
            # this in case you want to return different lags at once; don't care about that for now.
            # return np.correlate(x, y, mode="same")[data_size + 1: data_size + self.order + 1]  # many correlations at once
            # return np.correlate(x, y, mode="same")[data_size + lag] / (x.shape[0] - lag) - np.mean(x[lag:]) * np.mean(
            #     y[:-lag] if lag > 0 else y)  # x[:-0] does not work


class HakkarainenLorenzStatistics(LorenzLargerStatistics):
    """
    This class implements the statistics function from the Statistics protocol. This
    extracts the statistics following Hakkarainen et. al. [1] from the multivariate timesereis
    generated by solving Lorenz 95 odes.

    [1] J. Hakkarainen, A. Ilin, A. Solonen, M. Laine, H. Haario, J. Tamminen, E. Oja, and
    H. Järvinen. On closure parameter estimation in chaotic systems. Nonlinear Processes
    in Geophysics, 19(1):127–143, Feb. 2012.
    """

    def __init__(self, degree=2, cross=True):
        self.degree = degree
        self.cross = cross

        super(HakkarainenLorenzStatistics, self).__init__(moments=2, auto_covariances=1, covariances=1,
                                                          cross_covariances=1, degree=degree, cross=cross)
