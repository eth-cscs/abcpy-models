import re

import matplotlib.pyplot as plt
import numpy as np


def dist2(x, y):
    """Compute the square of the Euclidean distance between 2 arrays of same length"""
    return np.dot(x - y, x - y)


def contour_plot_no_grid(x, y, z, true_parameters, fig=None, ax=None, colorbar=False, log=False,
                         show_train_points=True, return_cntr=False, levels=None):
    """It produces a contour plot for data which is not on a grid. It performs some triangulation internally.
    Note that fig and ax have to be provided together."""

    if colorbar and ((fig is None) != (ax is None)):
        raise RuntimeError("Fig and Ax have to be provided together if you want to show colorbar!")

    if ax is None:
        ax = plt

    # apply log (if required)
    data = np.log(z) if log else z

    if levels is None:
        levels = 14

    # create plots
    ax.tricontour(x, y, data, levels=levels, linewidths=0.5, colors='k')
    cntr = ax.tricontourf(x, y, data, levels=levels, cmap="RdBu_r")  # , vmin=0, vmax=1)
    #    cntr = ax.tricontourf(x, y, data, levels=levels, cmap="RdYlGn")#, vmin=0, vmax=1)

    if colorbar:
        if fig is not None:
            fig.colorbar(cntr, ax=ax)
        else:
            ax.colorbar(cntr)

    if show_train_points:
        ax.plot(x, y, 'ko', ms=2)

    ax.plot(true_parameters[0], true_parameters[1], "ro", ms=5)
    if return_cntr:
        return cntr
    else:
        return


def estimate_params_from_distance(x, y, distances):
    return np.array([x[np.argmin(distances)], y[np.argmin(distances)]])


def read_simulations(file):
    f = open(file, 'r')
    param_line = True
    array_index = 0

    param1_list = []
    param2_list = []
    output_list = []  # these will contain the data
    simulation = np.zeros(72)  # each output has 72 points and 2 parameter values

    for x in f:
        # regex to extract the numbers
        line = re.findall(r"[-+]?[0-9]*[\.]?[0-9]+(?:[eE][-+]?[0-9]+)?", x)
        if param_line:
            param1_list.append(float(line[0]))
            param2_list.append(float(line[1]))
            param_line = False
            for i in range(2, len(line)):
                simulation[array_index] = float(line[i])
                array_index += 1
        else:
            for i in range(len(line)):
                simulation[array_index] = float(line[i])
                array_index += 1
                if array_index == 72:  # then we finished reading the array; reset and go to next simulation
                    output_list.append(simulation)
                    simulation = np.zeros(72)
                    param_line = True
                    array_index = 0

    param1 = np.array(param1_list)
    param2 = np.array(param2_list)
    output = np.array(output_list)

    return param1, param2, output


def compute_similarity_matrix(param1, param2, quantile=0.1, return_pairwise_distances=False):
    n_samples = len(param1)

    pairwise_distances = np.zeros([n_samples] * 2)

    for i in range(n_samples):
        for j in range(n_samples):
            pairwise_distances[i, j] = dist2(np.array([param1[i], param2[i]]), np.array([param1[j], param2[j]]))

    q10 = np.quantile(pairwise_distances[~np.eye(n_samples, dtype=bool)].reshape(-1), quantile)

    similarity_set = pairwise_distances < q10

    print("Fraction of similar pairs: ", np.sum(similarity_set) / n_samples ** 2)
    print("Fraction of similar pairs epurated by self-similarity: ",
          (np.sum(similarity_set) - n_samples) / n_samples ** 2)

    return (similarity_set, pairwise_distances) if return_pairwise_distances else similarity_set


# FUNCTIONS FOR ESTIMATION OF DIVERGENCES

def compute_distribution_from_distance(distances, beta=1, translate=True, normalize=True):
    """If translate is True, we translate the distances so that the minimum of them is 0, so that the exponential is 1.
    This may reduce numerical errors and have no impact, as the distribution is defined up to a multiplicative constant.
    If normalize is True, we then normalize the computed distribution on the samples, so that the sum over all samples
    is 1 (we essentialy get the P_i coefficients)"""

    if translate:
        distribution = np.exp(- beta * (distances - np.min(distances)))
    else:
        distribution = np.exp(- beta * distances)
    if normalize:
        distribution /= np.sum(distribution)

    return distribution


def estimate_norm_const(unnorm_distribution):
    return np.mean(unnorm_distribution)


def estimate_KL_distances(distances, true_distances, beta=1, tranlate=True):
    """We usually scale both distance to the same values before estimating the KL, otherwise the two distributions are
    differently concentrated."""
    # compute unnormalized distributions on the empirical samples
    distr = compute_distribution_from_distance(distances, beta=beta, translate=True, normalize=False)
    true_distr = compute_distribution_from_distance(true_distances, beta=beta, translate=True, normalize=False)
    # estimate norm constants:
    norm_const_true = estimate_norm_const(true_distr)
    norm_const = estimate_norm_const(distr)
    # compute KL estimate:
    return np.dot(beta * (true_distances - distances) + np.log(norm_const_true / norm_const), distr / norm_const)


def estimate_TV_distances(distances, true_distances, beta=1, tranlate=True):
    """We usually scale both distance to the same values before estimating the KL, otherwise the two distributions are
    differently concentrated."""
    # compute unnormalized distributions on the empirical samples
    distr = compute_distribution_from_distance(distances, beta=beta, translate=True, normalize=False)
    true_distr = compute_distribution_from_distance(true_distances, beta=beta, translate=True, normalize=False)
    # estimate norm constants:
    norm_const_true = estimate_norm_const(true_distr)
    norm_const = estimate_norm_const(distr)
    # compute KL estimate:
    return np.dot(np.abs((norm_const_true * distr) / (norm_const * true_distr) - 1), 0.5 * distr / norm_const)
