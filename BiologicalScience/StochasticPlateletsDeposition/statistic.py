from abcpy.statistics import Statistics

class Multiply(Statistics):
    """
    This class implements identity statistics not applying any transformation to the data, before the optional
    polynomial expansion step. If the data set contains n numpy.ndarray of length p, it returns therefore an
    nx(p+degree*p+cross*nchoosek(p,2)) matrix, where for each of the n points with p statistics, degree*p polynomial
    expansion term and cross*nchoosek(p,2) many cross-product terms are calculated.
    """

    def __init__(self, L, degree=1, cross=False, previous_statistics=None):
        """

        Parameters
        ----------
        degree : integer, optional
            Of polynomial expansion. The default value is 2 meaning second order polynomial expansion.
        cross : boolean, optional
            Defines whether to include the cross-product terms. The default value is True, meaning the cross product term
            is included.
        previous_statistics : Statistics class, optional
            It allows pipelining of Statistics. Specifically, if the final statistic to be used is determined by the
            composition of two Statistics, you can pass the first here; then, whenever the final statistic is needed, it
            is sufficient to call the `statistics` method of the second one, and that will automatically apply both
            transformations.
        """
        self.degree = degree
        self.cross = cross
        self.previous_statistics = previous_statistics
        self.L = L

    def statistics(self, data):
        """
        Parameters
        ----------
        data: python list
            Contains n data sets with length p.
        Returns
        -------
        numpy.ndarray
            nx(p+degree*p+cross*nchoosek(p,2)) matrix where for each of the n data points with length p,
            (p+degree*p+cross*nchoosek(p,2)) statistics are calculated.
        """

        # pipeline: first call the previous statistics:
        if self.previous_statistics is not None:
            data = self.previous_statistics.statistics(data)
        # the first of the statistics need to take list as input, in order to match the API. Then actually the
        # transformations work on np.arrays. In fact the first statistic transforms the list to array. Therefore, the
        # following code needs to be called only if the self statistic is the first, i.e. it does not have a
        # previous_statistic element.
        else:
            data = self._check_and_transform_input(data)
        # Chosen
        data = data[:, [6, 7, 8, 16, 17, 21, 22, 23]]
        # Expand the data with polynomial expansion
        result = self._polynomial_expansion(data)
        # Multiply with L
        result = result.dot(self.L.T)
        return result

class IdentityChosen(Statistics):
    """
    This class implements identity statistics not applying any transformation to the data, before the optional
    polynomial expansion step. If the data set contains n numpy.ndarray of length p, it returns therefore an
    nx(p+degree*p+cross*nchoosek(p,2)) matrix, where for each of the n points with p statistics, degree*p polynomial
    expansion term and cross*nchoosek(p,2) many cross-product terms are calculated.
    """

    def __init__(self, InformativeIndices = [3, 4, 6, 7, 8, 16, 17, 21, 22, 23], degree=1, cross=False, previous_statistics=None):
        """

        Parameters
        ----------
        degree : integer, optional
            Of polynomial expansion. The default value is 2 meaning second order polynomial expansion.
        cross : boolean, optional
            Defines whether to include the cross-product terms. The default value is True, meaning the cross product term
            is included.
        previous_statistics : Statistics class, optional
            It allows pipelining of Statistics. Specifically, if the final statistic to be used is determined by the
            composition of two Statistics, you can pass the first here; then, whenever the final statistic is needed, it
            is sufficient to call the `statistics` method of the second one, and that will automatically apply both
            transformations.
        """
        self.degree = degree
        self.cross = cross
        self.previous_statistics = previous_statistics
        self.InformativeIndices = InformativeIndices

    def statistics(self, data):
        """
        Parameters
        ----------
        data: python list
            Contains n data sets with length p.
        Returns
        -------
        numpy.ndarray
            nx(p+degree*p+cross*nchoosek(p,2)) matrix where for each of the n data points with length p,
            (p+degree*p+cross*nchoosek(p,2)) statistics are calculated.
        """

        # pipeline: first call the previous statistics:
        if self.previous_statistics is not None:
            data = self.previous_statistics.statistics(data)
        # the first of the statistics need to take list as input, in order to match the API. Then actually the
        # transformations work on np.arrays. In fact the first statistic transforms the list to array. Therefore, the
        # following code needs to be called only if the self statistic is the first, i.e. it does not have a
        # previous_statistic element.
        else:
            data = self._check_and_transform_input(data)
        # Chosen
        data = data[:, self.InformativeIndices]
        # Expand the data with polynomial expansion
        result = self._polynomial_expansion(data)
        return result