from abcpy.statistics import Statistics

class DiffusionIdentityStatistics(Statistics):
    """
    This class implements statistics for linear diffusion.
    """

    def __init__(self, degree=1, cross=False):
        self.degree = degree
        self.cross = cross

    def statistics(self, data):
        return data