# Distance learning for weakly supervised problems: 

In this folder you can find the files to run distance learning algorithms for the weakly supervised case of interest 
in the ABC setting we are considering.

There are some source files, defining functions and classes to run the algorithms, and some scripts to actually
run it. The volcanic model was run 400 times, and the resulting simulation-parameters pairs are provided in the file `training_data.txt`. 
The neural network based algorithms rely on the Pytorch package.
 The code for the contrastive and triplet loss algorithms has been adapted from https://github.com/adambielski/siamese-triplet.

Please check the `Volcano-distance-learning.ipynb` notebook to have more information on the different techniques, and to see
visualizations of the performance. 


## Source files: 

- `algorithms.py`: contains the training routines for the Mahalanobis distance learning algorithms SDML and the neural
 network-based ones with contrastive loss and triplet loss. It also contains the summary statistics learning approach based on neural 
 network embedding. 
- `datasets.py`: definitions of the Dataset objects for the neural network based algorithms.
- `losses.py`: definition of contrastive and triplet loss for distance learning using neural networks. 
- `networks.py`: contains the definition of the embedding networks used to transform the output of the volcanic inference
model, as well as the wrappers for running the contrastive and triplet loss algorithms (the `SiameseNet` is the wrapper to
be used with the contrastive loss, and the `TripletNet` is to be used with the triplet loss). It also contains a network version
used in the summary statistics learning approach, in which the size of the output of the network has to be the same as the
number of parameters you are inferring.
- `trainer.py`: contains the network training routine. 
- `utilities.py`: some utilities functions.

## Scripts: 

- `learn-networks.py`: runs the algorithms using neural networks on all of the training data and saves them
into the `../saved_networks` folder. These networks will be loaded by the ABC inference techniques and used to 
transform the output of the model. Also, the networks here are trained by defining the similarity sets with the best 
epsilon determined in the study of the epsilon sensitivity. 
- `cross_validation_script.py`: after splitting the 400 samples of the dataset into a training (300 samples) and test (100)
 dataset, it runs cross-validation in a Leave-One-Out fashion on the test set. Specifically, it first learn the distances 
 (or the summary statistics) on the
  training dataset, and then for each element in the test set, it is considered to be an observation once while all the 
  others serve as reference. Then, for each possible observation it computes the various distances from that to all elements
  in the reference set. It also stores the results in the folder `cross-valid-results`. As a very simple check, it computes also 
the parameters with the smallest distance value (amongst the reference ones) and the distances from the true parameters 
that originated the observation. This script uses epsilon corresponding to 10-th percentile. 
- `divergence_script.py`: starting from the results of the above, estimates the KL divergence and the Total Variation distance
between the probability distributions induced by each of the learned distances and the one induced by the true distance 
between the value of parameters, for each of the possible cross-validation splits on the test set. 
- `epsilon_sensitivity_NN_distances_script.py`: it performs the experiment about the sensitivity of the distance learning 
techniques to epsilon.  Specifically, it 
trains the network with contrastive and triplet loss for different choices of the threshold defining the pairwise similarity
set; then it estimates the KL divergence, in the same way as the above script using as observation each element in the test set, 
for each value of the threshold. Note that the statistics learning approach is not considered here as it does not depend on
the definition of a similarity set. The results are saved in the `cross-valid-results/epsilon_sensitivity_study` folder
- `visualization_distance_learning.py`: it produces the plots that can be found in the "Distance learning" section in the 
paper. In particular, it produces the distance contour for one possible cross-validation split of the test set and for
 each distance learning method considered in the paper; then, it produces histograms of the estimated KL divergence for 
 each of the methods, over the 100 possible Leave-One-Out splits of the test set. Finally, it draws the mean and standard deviation
 against the different values of the epsilon parameters used in the above script. The plots are saved in the `figures` folder. 
