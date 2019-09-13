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
#fake_data = volcano_model.forward_simulate([102.79663207, 98.75148678], 1, np.random.RandomState(), MPI.COMM_WORLD)
#print(fake_data[0].shape)

print("y_obs parameters : 173.87, 84.55")
obs_data = [
    np.array([15.18917688, 10.02628205,  5.73966907,  6.33214836,  9.90499409,
       6.33861612,  3.38936749,  5.40448485,  4.24736138,  4.24090759,
        7.92027537, 18.12765696,  3.73173343,  5.50238487,  5.18903677,
        6.10177494,  6.65568741, 10.33844081, 20.04962878, 21.9415658 ,
       19.84640171, 18.09792701, 20.74712285, 22.26741739, 20.99934227,
       19.62597728, 16.01909485, 17.14258782, 13.66039657, 16.17767254,
       13.7969041 , 13.23912777, 10.01644754, 10.01644754,  9.35800412,
       13.37838873,  9.17850549,  9.90499409,  9.90499409, 18.14852217,
        2.96388299, 20.090425  , 18.66639198, 19.0443569 , 15.83456207,
       12.78529858,  7.66271489,  3.81913262, 11.07588611, 22.14180809,
       26.43578794, 15.0581396 , 10.82687529,  2.92592224,  4.83855024,
       20.41073089, 19.71384789, 19.60073706, 12.56572833,  8.78751135,
       18.20711596, 16.94588622, 15.02197115,  9.46245009,  6.71842821,
        8.03805329, 11.9394705 ,  8.1754717 , 16.00248962, 17.49906809,
       16.72558858, 18.44861052]).reshape(-1, )]


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
#steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file = 6, 10, 1, 0.5, 0.03, 2.0, 1, None
journal = sampler.sample([obs_data], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file)
journal.save('VolcanojournalAPMCABC_simulated.jrnl')
