import numpy as np
from Model import Deposition
from Statistics import DepositionStatistics
from Distance import  DepositionDistance
from abcpy.distributions import Uniform, MultiNormal
from abcpy.inferences import SABC
from abcpy.output import Journal
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from abcpy.backends import BackendDummy
backend = BackendDummy()

#Define model 
lb, ub = [50,5,0.1,0.5e-3,0], [150,20,1.5,3e-3,10]
prior = Uniform(lb = [50,5,0.1,0.5e-3,0],ub = [150,20,1.5,3e-3,10])
pAd, pAg, pT, pF, aT = 110, 14.6, 0.6, 1.7e-3, 6
model = Deposition(prior, pAd, pAg, pT, pF, aT, seed = 1)
# Observed data
a =  np.array([ [0,0,0,172200,4808],
[20,1689,26.8,155100,1683],
[60,2004,29.9,149400,0],
[120,1968,31.3,140700,0],
[300,1946,36.6,125800,0]])

np.save('depo_experimental.npy',a)
BE = np.zeros(shape=(10,5))
index = 5
y_obs = np.load('depo_experimental.npy')[index]

# Define summary stat and distance
stat_calc = DepositionStatistics(degree = 1, cross = 0)
dist_calc = DepositionDistance(stat_calc)


mean = np.array([-13.0, .0, 7.0])
cov = np.eye(3)
kernel = MultiNormal(mean, cov, seed=1)

steps, epsilon, n_samples, n_samples_per_param = 2, 40, 1, 1
sampler = SABC(model, dist_calc, kernel, backend, seed = 1)
journal_sabc = sampler.sample([y_obs], steps, epsilon, n_samples, n_samples_per_param)
journal_sabc.save('experimental_5.jrnl')
samples = (journal_sabc.get_parameters(), journal_sabc.get_weights())
