# Spreading Process on Network

This folder contains models (Simple and Complex Contagion) and 
summary statistics, distance, Kernel, prior, infernce schemes 
to perform the Bayesian inference of the spreading processe parameters 
and the seed node. 

The subfolders SimpleContagion and ComplexContagion contains the needed 
functions for inference of the corresponding processes. 

The subfolder Networks contains 4 different networks: 1 Barabasi-Alber network
with 100 nodes, 1 Erdos-Renyi network with 100 nodes, Facebook social network
with 4039 nodes and Indian village contact network with 354 nodes.

The SimpleContagionExample.py and ComplexContagionExample.py illustrates how 
to simulate spreading processes on a chosen network and how to infer parameters
using that simulated epidemic. 

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

