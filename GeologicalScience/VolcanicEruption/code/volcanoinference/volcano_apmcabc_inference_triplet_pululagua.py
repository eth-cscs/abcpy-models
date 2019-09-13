import numpy as np
from fileio import *
from mpi4py import MPI
import argparse
from abcpy.backends import BackendMPI as Backend

import logging
logging.basicConfig(level=logging.INFO)

print('Reading process per model')

parser = argparse.ArgumentParser()
parser.add_argument("npm", help="number of mpi process per model", type=int)
args = parser.parse_args()
process_per_model = args.npm
if process_per_model >=MPI.COMM_WORLD.Get_size():
    raise "number of process per model must be lower than number of MPI process (one process is dedicated to the scheduler)"
backend = Backend(process_per_model=process_per_model)

print('backend Initiated')

from abcpy.continuousmodels import Uniform
from Model import Volcano
U0 = Uniform([[100], [300]], 'U0')
L = Uniform([[30], [100]], 'L')

# define the model
volcano_model = Volcano([U0, L], name='volcano_model')
print("y_obs from Pululagua eruption")
obs_data = [
 np.array([70.0, 100.0, 170.0, 150.0, 290.0, 220.0, 235.0, 180.0, 275.0, 175.0, 170.0,
165.0, 140.0, 100.0, 100.0, 120.0, 110.0, 80.0, 25.0, 20.0, 30.0, 40.0, 30.0,
20.0, 15.0, 30.0, 70.0, 70.0, 100.0, 70.0, 75.0, 75.0, 100.0, 110.0, 115.0, 80.0,
125.0, 280.0, 350.0, 350.0, 180.0, 25.0, 40.0, 40.0, 55.0, 70.0, 175.0, 260.0,
200.0, 140.0, 50.0, 10.0, 120.0, 300.0, 250.0, 60.0, 90.0, 105.0, 150.0, 170.0,
120.0, 110.0, 120.0, 190.0, 190.0, 130.0, 130.0, 110.0, 70.0, 90.0, 85.0, 65.0]).reshape(-1, )]

# define statistics
# from abcpy.statistics import Identity
# statistics_calculator = Identity(degree=1, cross=False)

from Statistics import NeuralEmbeddingStatistics
from Statistics import load_net
from distance_learning.networks import EmbeddingNet

# Here we load the network and pass it to the statistics calculator
embedding_net_triplet = load_net("saved-networks/triplet.pth", EmbeddingNet)
embedding_net_triplet.eval()
statistics_calculator = NeuralEmbeddingStatistics(embedding_net_triplet, degree=1, cross=False)


# define distance
from Distance import DistanceVolcano
distance_calculator = DistanceVolcano(statistics_calculator)

#print(distance_calculator.distance(fake_data, obs_data))

# # define sampling scheme
# # (seed is not used for now)
from abcpy.inferences import APMCABC
sampler = APMCABC([volcano_model], [distance_calculator], backend, seed=1)
print('sampling')
steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file = 6, 100, 1, 0.1, 0.03, 2.0, 1, None
#steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file = 1, 100, 1, 0.1, 0.03, 2.0, 1, 'VolcanojournalAPMCABC_pululagua.jrnl'
journal = sampler.sample([obs_data], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file)
journal.save('VolcanojournalAPMCABC_pululagua.jrnl')
