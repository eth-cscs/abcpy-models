import numpy as np
from abcpy.distributions import Distribution

class DiffusionPrior(Distribution):
    """
    This class implements prior distribution for linear diffusion.
    """
    """
    This class implements a p-dimensional uniform Prior distribution in a closed interval.
    """

    def __init__(self, index, network, vec, lb: np.ndarray, ub: np.ndarray, seed=None):
        """
        Defines the upper and lower bounds of a p-dimensional uniform Prior distribution in a closed interval.

        Parameters
        ----------
        index: boolean
            Indicating whether discrete or continuous uniform, if 1 then continous uniform, ow discrete uniform.
        vec
        lb: numpy.ndarray or a list
            Vector containing p lower bounds
        ub: numpy.ndarray or a list
            Vector containing p upper bounds
        seed: integer
            Initial seed for the random number generator

        """
        self.vec = vec
        self.lb, self.ub = self._check_parameters(lb, ub)
        self.index = index
        self.rng = np.random.RandomState(seed)
        # self.weight = np.zeros(shape=(len(vec),))
        # for ind in range(len(vec)):
        #    self.weight[ind] = 1/network.degree(vec[ind])
        # self.weight = self.weight/sum(self.weight)

    def set_parameters(self, params):
        lb = params[0]
        ub = params[1]
        self.lb, self.ub = self._check_parameters(lb, ub)

    def reseed(self, seed):
        self.rng.seed(seed)

    def sample(self, k):
        samples = np.zeros(shape=(k, len(self.lb)))
        for j in range(0, len(self.lb)):
            if self.index[j] == 1:
                samples[:, j] = self.rng.uniform(self.lb[j], self.ub[j], k)
            else:
                samples[:, j] = self.vec[self.rng.randint(0, len(self.vec), k)]
                # samples[:,j] = self.rng.choice(self.vec,size = k, p = self.weight)
        return samples

    def pdf(self, x):
        pdf_value = 1
        for j in range(0, len(self.lb)):
            if self.index[j] == 1:
                pdf_value = pdf_value * (1 / (self.ub[j] - self.lb[j])) \
                            * np.product(np.greater_equal(x[j], self.lb[j]) * np.less_equal(x[j], self.ub[j]))
            else:
                pdf_value = pdf_value * (1 / len(self.vec)) * (len(np.where(self.vec == x[j])[0]) > 0)
                # pdf_value = pdf_value*(self.weight[np.where(self.vec==x[j])[0]])*(len(np.where(self.vec==x[j])[0])>0)
        return pdf_value

    def _check_parameters(self, lb, ub):
        new_lb = new_ub = None
        if isinstance(lb, (list, np.ndarray)):
            new_lb = np.array(lb)
        else:
            raise TypeError('The lower bound is not of allowed types')

        if isinstance(ub, (list, np.ndarray)):
            new_ub = np.array(ub)
        else:
            raise TypeError('The upper bound is not of allowed types')

        if new_lb.shape != new_ub.shape:
            raise BaseException('Dimension of lower bound and upper bound is not same.')

        return (new_lb, new_ub)