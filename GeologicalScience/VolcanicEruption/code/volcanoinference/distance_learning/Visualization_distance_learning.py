# add to the pythonpath the current folder (otherwise it does not find the script files)
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir)))

import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import train_test_split
from utilities import read_simulations, contour_plot_no_grid

# load distance files

best_epsilon = 0.6

embedding_distances_true = np.load("cross-valid-results/true_embedding_distances.npy")
embedding_distances_naive = np.load("cross-valid-results/naive_embedding_distances.npy")
embedding_distances_sdml = np.load("cross-valid-results/sdml_embedding_distances.npy")
embedding_distances_contrastive = np.load("cross-valid-results/contrastive_batch32_epochs400_embedding_distances.npy")
embedding_distances_triplet = np.load("cross-valid-results/triplet_batch16_epochs800_embedding_distances.npy")
embedding_distances_FP_nn = np.load("cross-valid-results/FP_nn_epochs400_embedding_distances.npy")
embedding_distances_triplet_best_eps = np.load("cross-valid-results/epsilon_sensitivity_study/"
                                               "triplet_batch16_epochs800_embedding_distances"
                                               + "epsilon_{:.2f}".format(best_epsilon) + ".npy")

# define output folder
folder = "figures/"

# load parameter: 
param1, param2, output = read_simulations('training_data.txt')
output_train, test_data, train_labels, test_labels = train_test_split(output, np.column_stack((param1, param2)),
                                                                      test_size=100, shuffle=False)

# fix the index of the simulation you want to plot and filter the data:

index = 65
param1_obs = test_labels[index, 0]
param2_obs = test_labels[index, 1]
true_parameters = np.array([param1_obs, param2_obs])

param1_ref = test_labels[np.arange(len(test_labels)) != index, 0]
param2_ref = test_labels[np.arange(len(test_labels)) != index, 1]

distances_true = embedding_distances_true[index]
distances_naive = embedding_distances_naive[index]
distances_sdml = embedding_distances_sdml[index]
distances_contrastive = embedding_distances_contrastive[index]
distances_triplet = embedding_distances_triplet[index]
distances_FP_nn = embedding_distances_FP_nn[index]
distances_triplet_best_eps = embedding_distances_triplet_best_eps[index]


# PLOT THE DISTANCE CONTOURS:
def rescale(distances):
    return distances / np.max(distances)


log = False
show_train_points = True
levels = 20

title_size = 15
size = 18

titles = ["Euclidean distances", "SDML distances", "Contrastive loss distances", "Triplet loss distances",
          "True distances", "Semiautomatic NN distances", r"Triplet loss distances (best $\epsilon$)"]
results = [distances_naive, distances_sdml, distances_contrastive, distances_triplet, distances_true, distances_FP_nn,
           distances_triplet_best_eps]

for title, result in zip(titles, results):
    fig, ax = plt.subplots(figsize=(16 / 5, 4))

    ax.set_title(title, size=title_size)

    ax.set_ylabel(r"$U_0\ [m/s]$", size=size)

    ax.set_xlabel(r"$R_0\ [m]$", size=size)

    # fig.tight_layout()
    
    colorbar = True if title in ["Semiautomatic NN distances", r"Triplet loss distances (best $\epsilon$)"] else False

    cntr = contour_plot_no_grid(param2_ref, param1_ref, rescale(np.transpose(result)), true_parameters[::-1], ax=ax,
                                log=log, show_train_points=show_train_points, fig=fig, levels=levels, return_cntr=True,
                                colorbar=colorbar)

    plt.savefig(folder + "distance_maps_" + title.replace(' ', '_') + '.eps', bbox_inches='tight')

# save now colorbar separately: 

fig, ax = plt.subplots(figsize=(0.35, 4))
fig.colorbar(cntr, cax=ax)

plt.savefig(folder + "distance_maps_colorbar.eps", bbox_inches='tight')


# PLOT THE DIVERGENCE HISTOGRAMS, IN CROSS-VALIDATION SETTING:
def write_mean_std(ax, vector):
    # exclude the nan elements
    ax.text(0.5, 0.9, r"$m$={:.2f}".format(np.median(vector[np.isfinite(vector)])), size=20, transform=ax.transAxes)
    ax.text(0.5, 0.8, r"$\mu$={:.2f}".format(np.mean(vector[np.isfinite(vector)])), size=20, transform=ax.transAxes)
    ax.text(0.5, 0.7, r"$\sigma$={:.2f}".format(np.std(vector[np.isfinite(vector)])), size=20, transform=ax.transAxes)


def plot_hist(ax, vector, alpha=0.7, rwidth=0.85, alpha_grid=0.75, bins=None):
    ax.grid(axis='y', alpha=alpha_grid)
    ax.hist(np.array(vector), density=True, color='#0504aa',
            alpha=alpha, rwidth=rwidth, bins=bins)


def compute_distance_pair_params(vector):
    return np.sqrt(np.sum(np.power(vector, 2), axis=1))


# load data
KL_IS_naive = np.load("cross-valid-results/naive_estimated_KL_IS_distance_functions.npy")
KL_IS_sdml = np.load("cross-valid-results/sdml_estimated_KL_IS_distance_functions.npy")
KL_IS_contrastive = np.load("cross-valid-results/contrastive_batch32_epochs400_estimated_KL_IS_distance_functions.npy")
KL_IS_triplet = np.load("cross-valid-results/triplet_batch16_epochs800_estimated_KL_IS_distance_functions.npy")
KL_IS_FP_nn = np.load("cross-valid-results/FP_nn_epochs400_estimated_KL_IS_distance_functions.npy")
KL_IS_triplet_best_eps = np.load("cross-valid-results/epsilon_sensitivity_study/"
                                 "triplet_batch16_epochs800_estimated_KL_IS_distance_functions"
                                 + "epsilon_{:.2f}".format(best_epsilon) + ".npy")

bins = np.linspace(0, 11, 20)

title_size = 16
size = 14

titles = ["Euclidean distances", "SDML distances", "Contrastive loss distances", "Triplet loss distances",
          "Semiautomatic NN distances", r"Triplet loss distances (best $\epsilon$)"]
results = [KL_IS_naive, KL_IS_sdml, KL_IS_contrastive, KL_IS_triplet, KL_IS_FP_nn, KL_IS_triplet_best_eps]

for title, result in zip(titles, results):
    fig, ax = plt.subplots(figsize=(4 * 0.9, 4 * 0.9))
    ax.set_xlim([-1, 12])
    ax.set_ylim([0, .8])

    ax.set_title(title, size=title_size)

    ax.set_ylabel("Density", size=size)

    ax.set_xlabel("Estimated KL divergence", size=size)

    fig.tight_layout()

    write_mean_std(ax, result)
    plot_hist(ax, result, bins=bins)

    plt.savefig(folder + "divergence-histogram_" + title.replace(' ', '_') + '.eps', bbox_inches="tight")

# PLOT THE SENSITIVITY OF THE CONTRASTIVE AND TRIPLET TECHNIQUE TO THE CHOICE OF EPSILON, ACCORDING TO THE ESTIMATED KL DIVERGENCE
epsilon_values = np.arange(0.02, 1.02, 0.02)
folder_eps_data = "cross-valid-results/epsilon_sensitivity_study"

# load the files:

KL_list_contr = []
KL_list_tripl = []
KL_list_sdml = []

epsilon_accepted_contr = []
epsilon_accepted_tripl = []
epsilon_accepted_sdml = []

for epsilon in epsilon_values:
    try:
        KL_divergence_contr = np.load(folder_eps_data + "/" + 'contrastive_batch32_epochs400' + "_estimated_KL_IS_distance_functions" + 
                                      "epsilon_{:.2f}".format(epsilon) + ".npy")
        KL_list_contr.append(KL_divergence_contr)
        epsilon_accepted_contr.append(epsilon)
    except FileNotFoundError:
        print("Epsilon {:.2f} not working for contrastive".format(epsilon))
        
    try:
        KL_divergence_tripl = np.load(folder_eps_data + "/" + 'triplet_batch16_epochs800' + "_estimated_KL_IS_distance_functions" + 
                                      "epsilon_{:.2f}".format(epsilon) + ".npy")
        KL_list_tripl.append(KL_divergence_tripl)
        epsilon_accepted_tripl.append(epsilon)
    except FileNotFoundError:
        print("Epsilon {:.2f} not working for triplet".format(epsilon))
    
    try:
        KL_divergence_sdml = np.load(folder_eps_data + "/" + 'sdml' + "_estimated_KL_IS_distance_functions" + 
                                     "epsilon_{:.2f}".format(epsilon) + ".npy")
        KL_list_sdml.append(KL_divergence_sdml)
        epsilon_accepted_sdml.append(epsilon)
    except FileNotFoundError:
        print("Epsilon {:.2f} not working for SDML".format(epsilon))
        
KL_list_contr = np.array(KL_list_contr).T
KL_list_tripl = np.array(KL_list_tripl).T
KL_list_sdml = np.array(KL_list_sdml).T

def set_colors(bp, color1, color2, alpha):

    for box in bp['boxes']:
        # change outline color
        box.set( color=color1, alpha=alpha)
        # change fill color
        #box.set( facecolor = '#1b9e77' )

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color=color1, alpha=alpha)#, linewidth=2)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color=color1, linewidth=2, alpha=alpha)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color=color2, linewidth=2)
    
    ## change color and linewidth of the medians
    for mean in bp['means']:
        mean.set(color=color2, linewidth=2)

        
    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='.', alpha=alpha, markerfacecolor=color1, markeredgecolor=color1)
        
alpha = 0.5
contr_color1 = 'tomato'
contr_color2 = 'red' 
tripl_color1 = 'royalblue'
tripl_color2 = 'blue' 
sdml_color1 = 'limegreen'
sdml_color2 = 'green' 
mean = False

fig, ax = plt.subplots(figsize=(9,3))

bp_contr = ax.boxplot(KL_list_contr, notch=False, showmeans=mean, meanline=True, showfliers=True, patch_artist=True,
                positions=epsilon_accepted_contr, widths = 0.015, )
set_colors(bp_contr, contr_color1, contr_color2, alpha)
    
bp_tripl = ax.boxplot(KL_list_tripl, notch=True, showmeans=mean, meanline=True, showfliers=True, patch_artist=True,
                  positions=epsilon_accepted_tripl, widths = 0.015, )
set_colors(bp_tripl, tripl_color1, tripl_color2, alpha)

bp_sdml = ax.boxplot(KL_list_sdml, notch=True, showmeans=mean, meanline=True, showfliers=True, patch_artist=True,
                  positions=epsilon_accepted_sdml, widths = 0.015, )
set_colors(bp_sdml, sdml_color1, sdml_color2, alpha)

ax.legend([bp_sdml["boxes"][0], bp_contr["boxes"][0], bp_tripl["boxes"][0]],
          ['SDML', 'Contrastive', 'Triplet'])

ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ax.set_xticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
ax.set_xticks(epsilon_values, minor = True)
ax.set_xlim([-0.03,1.05])

ax.set_ylabel("Estimated KL divergence")
ax.set_xlabel("Quantile")

fig.savefig('figures/epsilon_sensitivity_NN_boxplot.pdf', bbox_inches="tight")
