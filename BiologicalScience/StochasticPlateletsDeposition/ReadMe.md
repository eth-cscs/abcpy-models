# Platelets Deposition Model
Numerical model of platelet deposition and adhesion in human blood as described in Chopard (2017). 
The summary statistics, distance functions and the ABC inference schemes are the same ones used in 
Dutta et. al. (2018).

## Reference
- Chopard B, de Sousa DR, Lätt J, Mountrakis L, Dubois F, Yourassowsky C, et al. A physical
description of the adhesion and aggregation of platelets. Royal Society Open Science 4 (2017) 170219.
- Dutta R, Chopard B, Lätt J, Dubois F, Zouaoui Boudjeltia K, Mira, A Parameter estimation of 
platelets deposition: Approximate Bayesian computation with high performance computing. 
Frontiers in Physiology: Computational Physiology and Medicine (2018)

## Dependencies 
- g++ 
- Swig 
- abcpy0.5.2

## Setup 

- Run make 
- You can run an example ABC inference scheme on the platelet deposition model by 'python3 Inference.py'
- (if MPI is installed) 
1. Modify line 31 of Inference.py from 'from abcpy.backends import BackendDummy as Backend' to 'from abcpy.backends import BackendMPI as Backend'
2. To run: 'mpirun -np 8 python3 Inference.py'
