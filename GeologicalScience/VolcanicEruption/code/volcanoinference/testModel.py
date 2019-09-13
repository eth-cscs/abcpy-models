import numpy as np
from fileio import *
from mpi4py import MPI
import argparse
from abcpy.backends import BackendMPI as Backend

import logging
logging.basicConfig(level=logging.INFO)

from abcpy.continuousmodels import Uniform
from Model import Volcano
U0 = Uniform([[100], [300]], 'U0')
L = Uniform([[30], [100]], 'L')

# define the model
volcano_model = Volcano([U0, L], name='volcano_model')
fake_data = volcano_model.forward_simulate([102.79663207, 98.75148678], 1, np.random.RandomState(), MPI.COMM_WORLD)
print(fake_data[0].shape)
print(fake_data)


