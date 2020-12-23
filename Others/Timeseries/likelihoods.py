import numpy as np
from scipy.stats import multivariate_normal


def ma2_likelihood(x, theta_1, theta_2):
    cov_matrix = np.zeros(x.shape * 2)
    a = 1 + theta_1 ** 2 + theta_2 ** 2
    cov_matrix[0, 0] = 1
    cov_matrix[1, 1] = 1 + theta_1 ** 2
    b = theta_1 * (1 + theta_2)
    cov_matrix[0, 1] = theta_1
    cov_matrix[1, 0] = theta_1
    for i in range(2, x.shape[0]):
        cov_matrix[i, i] = a
        cov_matrix[i - 1, i] = b
        cov_matrix[i, i - 1] = b
        cov_matrix[i - 2, i] = theta_2
        cov_matrix[i, i - 2] = theta_2
    return multivariate_normal.logpdf(x, cov=cov_matrix)


def ar2_likelihood2(x, theta_1, theta_2):
    # as above but faster
    means = np.zeros_like(x)
    means[1] = theta_1 * x[0]
    means[2:] = theta_1 * x[1:-1] + theta_2 * x[0:-2]

    return multivariate_normal.logpdf(x, mean=means)
