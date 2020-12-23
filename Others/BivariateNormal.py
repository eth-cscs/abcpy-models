import numpy as np
from abcpy.continuousmodels import ProbabilisticModel, Continuous, InputConnector
from scipy.stats import multivariate_normal


# following definition in https://mathworld.wolfram.com/BivariateNormalDistribution.html. We also give the loglikelihood
# and the unnormalized loglikelihood. Notice that this model produces n_samples iid samples and flattens them out.

class BivariateNormal(ProbabilisticModel, Continuous):

    def __init__(self, parameters, n_samples=10, name='Bivariate_Normal'):
        input_parameters = InputConnector.from_list(parameters)
        self.n_samples = n_samples
        super(BivariateNormal, self).__init__(input_parameters, name)

    @staticmethod
    def _create_mean_cov(mu1, mu2, sigma1, sigma2, rho):
        mean = np.array([mu1, mu2])
        off_diag_elem = rho * sigma1 * sigma2
        cov = np.array([[sigma1 ** 2, off_diag_elem], [off_diag_elem, sigma2 ** 2]])
        return mean, cov

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        result = self.sample(*input_values, rng=rng)
        return [result]

    def sample(self, mu1, mu2, sigma1, sigma2, rho, rng=None):
        mean, cov = self._create_mean_cov(mu1, mu2, sigma1, sigma2, rho)
        x = multivariate_normal.rvs(mean=mean, cov=cov, size=self.n_samples, random_state=rng)
        return x.flatten()

    def logpdf(self, x, mu1, mu2, sigma1, sigma2, rho):
        x = x.reshape(self.n_samples, 2)
        mean, cov = self._create_mean_cov(mu1, mu2, sigma1, sigma2, rho)
        logpdfs = multivariate_normal.logpdf(x=x, mean=mean, cov=cov, allow_singular=True)
        return np.sum(logpdfs)

    def logpdf_unnorm(self, x, mu1, mu2, sigma1, sigma2, rho):
        logpdfs = self.logpdf(x, mu1, mu2, sigma1, sigma2, rho)
        logpdfs += np.log(sigma1) + np.log(sigma2) + 0.5 * np.log(1 - rho ** 2)
        return np.sum(logpdfs)

    def get_output_dimension(self):
        return 2 * self.n_samples

    def _check_input(self, values):
        return True

    def _check_output(self, values):
        return True
