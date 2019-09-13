import numpy as np
import torch
from sklearn.model_selection import LeaveOneOut, train_test_split
from sklearn.preprocessing import MinMaxScaler
from scipy.interpolate import Rbf
from tqdm import tqdm

from utilities import read_simulations, estimate_TV_distances, estimate_KL_distances

param1, param2, output = read_simulations('training_data.txt')
# instantiate scaler
scaler = MinMaxScaler()

# should insert loop over all possible techniques for which we have the embedding distances.
folder = "cross-valid-results"
# we exclude now sdml because we were not able to fit all the possible splits and we lost the information on which true
# parameter value was successful.
names = ("true", "naive", "sdml", "contrastive_batch32_epochs400", "triplet_batch16_epochs800", "FP_nn_epochs400")
names = ["FP_nn_epochs400"]

output_train, test_data, train_labels, test_labels = train_test_split(output, np.column_stack((param1, param2)),
                                                                      test_size=100, shuffle=False)
param1_train = train_labels[:, 0]
param2_train = train_labels[:, 1]
param1_test = test_labels[:, 0]
param2_test = test_labels[:, 1]

# open the ground truth (i.e. the tre distances on those points)
embedding_distances_true = np.load(folder + "/true_embedding_distances.npy")

for name in names:

    # flag for sdml case
    is_sdml = True if name == "sdml" else False

    print(name)

    loo = LeaveOneOut()
    estimated_KL = []
    estimated_TV = []

    # load the embedding distances
    embedding_distances = np.load(folder + "/" + name + "_embedding_distances.npy")

    if is_sdml:
        true_parameters_successful = np.load(folder + "/" + name + "_true_params.npy")
        sdml_index = 0

    for ref_index, obs_index in tqdm(loo.split(test_data, test_labels)):  # this returns the indices

        # split the dataset
        n_samples = len(ref_index)

        output_ref = test_data[ref_index]
        param1_ref = param1_test[ref_index]
        param2_ref = param2_test[ref_index]

        observation = test_data[obs_index]
        param1_obs = param1_test[obs_index]
        param2_obs = param2_test[obs_index]
        true_parameters = np.array([param1_obs, param2_obs]).reshape(-1)

        if is_sdml:
            if param1_obs == true_parameters_successful[sdml_index, 0] and param2_obs == true_parameters_successful[
                sdml_index, 1]:
                sdml_index += 1
            else:
                continue  # skip this split as it was not successful when fitting the sdml technique. 

        test_index = sdml_index - 1 if is_sdml else obs_index[0]

        # select the correct distance now: 
        distances = embedding_distances[test_index]
        true_distances = embedding_distances_true[test_index]

        # scale them in [0,1]
        distances = scaler.fit_transform(distances.reshape(-1, 1)).reshape(-1)
        true_distances = scaler.fit_transform(true_distances.reshape(-1, 1)).reshape(-1)

        # estimate the KL and append
        estimated_KL.append(estimate_KL_distances(distances, true_distances, beta=1))
        estimated_TV.append(estimate_TV_distances(distances, true_distances, beta=1))

    np.save(folder + "/" + name + "_estimated_KL_IS_distance_functions", np.array(estimated_KL))
    np.save(folder + "/" + name + "_estimated_TV_IS_distance_functions", np.array(estimated_TV))
