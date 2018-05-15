# Spreading Process on Network

This folder contains models (Simple and Complex Contagion) and 
summary statistics, distance, Kernel, prior, infernce schemes 
to perform the Bayesian inference of the spreading processe parameters 
and the seed node. Further we also include some codes to compute Bayes Estimate, 
visualize posterior distribution etc. 

The subfolders SimpleContagion and ComplexContagion contains the needed 
functions for inference of the corresponding processes. 

The subfolder Networks contains 4 different networks: Simulated Barabasi-Alber network
with 100 nodes, Simulated Erdos-Renyi network with 100 nodes, Facebook social network
with 4039 nodes and Indian village contact network with 354 nodes.

# Folders and Files:
--------------------------------
* SimpleContagion: This subfolder contains the needed functions for inference of the corresponding processes. 
* ComplexContagion: This subfolder contains the needed functions for inference of the corresponding processes. 
* SimpleContagionExample.py and ComplexContagionExample.py: Simulates corresponding spreading processes on a 
chosen network and infers parameters for that simulated epidemic. 
* BayeEstimate.py: Computes the Bayes Estimates
* SimpleContagionBEVisualization.py and ComplexContagionBEVisualization.py: Compute the Bayes Estimates for 
the 100 simulated epidemics on the simulated networks and creates the corresponding Figures.
* SimpleContagionPosteriorVisualization.py and ComplexContagionPosteriorVisualization.py: Creates the figure 
of the posterior distributions.
* Results: This folder is used to save the simulated epidemics and inferred results. (Empty, contains are avaiable on request)
* Figure: This folder is used to save the Figures generated from the visualization file. (Empty, contains are avaiable on request)
 

## Reference
The above models and inference schemes are described in detail in the manuscript: 
R. Dutta, J. P. Onnela, A. Mira, "Bayesian Inference of Spreading Processes on 
Networks", 2017, arXiv:1709.08862


## Dependencies 
- Networkx Python package 
- abcpy0.3.0

## Setup 
- tar xzvf Networks.tar.gz
- python3 SimpleContagionExample.py or python3 ComplexContagionExample.py

