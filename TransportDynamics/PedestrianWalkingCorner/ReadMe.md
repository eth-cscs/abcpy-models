# Pedestrian walking through a corner
This floor field model is for simulating pedestrians' walking throug a corner.
Three models are taken to generate the static floor field which is responsible to governing pedestrian navigation:
- a simple approach proposed by Katsuhiro Nishinari
- the block algorithm proposed by Li et al.
- the function proposed by Dias and Lovreglio

## Reference
1. Li, Shengnan, et al. "Block-based floor field model for pedestrian’s walking through corner." Physica A: Statistical Mechanics and its Applications 432 (2015): 337-353.
2. Dias, Charitha and Lovreglio, Ruggiero "Calibrating cellular automaton models for pedestrians walking through corners." Physics Letters A 382.19 (2018): 1255-1261.

## Dependencies 
- abcpy0.5.1

## Setup
You can run an example ABC inference scheme on this model by running 'python3 Inference.py'
