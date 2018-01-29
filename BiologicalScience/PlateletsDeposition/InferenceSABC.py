import numpy as np
from examples.bastien.param_5.Deposition_model import Deposition, DepositionStatisticsCombined, DepositionDistanceCombined
#DepositionRelativeErrorDistance, DepositionIdentityStatistics, DepositionDistance, DepositionStatistics
from abcpy.distributions import Uniform, MultiNormal
from abcpy.inferences import SABC
from abcpy.output import Journal
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from abcpy.backends import BackendDummy
backend = BackendDummy()

#==============================================================================
#import findspark
#findspark.init()
#import pyspark
#sc = pyspark.SparkContext()
#from abcpy.backends import BackendSpark as Backend
#backend = Backend(sc, parallelism=4)
# 
#from abcpy.backend_mpi import BackendMPI as Backend
#backend = Backend()
# 
#==============================================================================
#Define model 
lb, ub = [50,5,0.1,0.5e-3,0], [150,20,1.5,3e-3,10]
prior = Uniform(lb = [50,5,0.1,0.5e-3,0],ub = [150,20,1.5,3e-3,10])
pAd, pAg, pT, pF, aT = 110, 14.6, 0.6, 1.7e-3, 6
model = Deposition(prior, pAd, pAg, pT, pF, aT, seed = 1)
# Simulate observed data
#==============================================================================
#==============================================================================
#import time
#start_time = time.time()
#y_obs = model.simulate(1)
#print("--- %s seconds ---" % (time.time() - start_time))
# # np.save('trial_data.npy',y_obs)
#==============================================================================

#a =  np.array([ [0,0,0,172200,4808],
#[20,1689,26.8,155100,1683],
#[60,2004,29.9,149400,0],
#[120,1968,31.3,140700,0],
#[300,1946,36.6,125800,0]])

#np.save('depo_experimental.npy',a)
# 
#==============================================================================
#y_obs = np.load('trial_data.npy')
BE = np.zeros(shape=(10,5))
index = 5
y_obs = np.load('depo_experimental.npy')[index]

# Define summary stat and distance
stat_calc = DepositionStatisticsCombined(degree = 1, cross = 0)
dist_calc = DepositionDistanceCombined(stat_calc)

#stat_calc = DepositionStatistics(degree = 1, cross = 0)
#dist_calc = DepositionDistance(stat_calc)
#stat_calc = DepositionIdentityStatistics(degree = 1, cross = 0)
#dist_calc = DepositionRelativeErrorDistance(stat_calc)


mean = np.array([-13.0, .0, 7.0])
cov = np.eye(3)
kernel = MultiNormal(mean, cov, seed=1)

# use the rejection sampling scheme for Inference 
#sampler = RejectionABC(model, dist_calc, backend, seed = 1)
#n_samples, n_samples_per_param, epsilon = 4, 1, 10
#journal = sampler.sample(y_obs, n_samples, n_samples_per_param, epsilon)
#samples=journal.get_parameters()

#import time
#start_time = time.time()
#steps, epsilon, n_samples, n_samples_per_param = 2, 40, 1, 1
#sampler = SABC(model, dist_calc, kernel, backend, seed = 1)
#journal_sabc = sampler.sample([y_obs], steps, epsilon, n_samples, n_samples_per_param)
#journal_sabc.save('examples/bastien/Result/experimental_5.jrnl')
#samples = (journal_sabc.get_parameters(), journal_sabc.get_weights())
#print("--- %s seconds ---" % (time.time() - start_time))

from scipy.integrate import simps
import pymc3
jj = Journal.fromFile('examples/bastien/Result/experimental_5_'+str(index)+'.jrnl')
samples = jj.get_parameters()
CI = pymc3.stats.hpd(samples)
label1 = [r'$p(p_{Ad}|x): empirical$', r'$p(p_{Ag}|x): empirical$', r'$p(p_{T}|x): empirical$', r'$p(p_{F}|x): empirical$', r'$p(a_{T}|x): empirical$']
label1s = [r'$p(p_{Ad}|x): Smoothed \ posterior \ density$', r'$p(p_{Ag}|x): Smoothed \ posterior \ density$', r'$p(p_{T}|x): Smoothed \ posterior \ density$', r'$p(p_{F}|x): Smoothed \ posterior \ density$', r'$p(a_{T}|x): Smoothed \ posterior \ density$']
label2 = [r'$p_{Ad}^{man. est.}$',r'$p_{Ag}^{man. est.}$',r'$p_{T}^{man. est.}$',r'$p_{F}^{man. est.}$',r'$a_{T}^{man. est.}$']
label3 = [r'$E(p_{Ad}|x)$: Bayes \ estimate', r'$E(p_{Ag}|x): Bayes \ estimate$', r'$E(p_{T}|x): Bayes \ estimate$', r'$E(p_{F}|x): Bayes \ estimate$', r'$E(a_{T}|x): Bayes \ estimate$']
label4 = [r'$p_{Ad}$',r'$p_{Ag}$', r'$p_{T}$', r'$p_{F}$' , r'$a_{T}$']
theta = [110, 14.6, 0.6, 1.7e-3,6]
Bayes_estimate = np.mean(samples, axis=0)
BE[index,:] = Bayes_estimate
for ind in range(5): 
    #xmin, xmax = max(samples[:,ind]), min(samples[:,ind])
    xmin, xmax = lb[ind], ub[ind]
    positions = np.linspace(xmin, xmax, samples.shape[0])
    gaussian_kernel = gaussian_kde(samples[:,ind].reshape(samples.shape[0],),bw_method=0.01)
    values = gaussian_kernel(positions)
    gaussian_kernel = gaussian_kde(samples[:,ind].reshape(samples.shape[0],),bw_method=.1)
    values1 = gaussian_kernel(positions)
    #Bayes_estimate[ind] = positions[list(values1).index(max(values1))]
    #simps(positions*values1,positions,axis=-1)
    plt.figure()
    #plt.plot(positions,values, label = label1[ind])
    plt.plot(positions,values1, label = label1s[ind], color = 'b')
    plt.plot([Bayes_estimate[ind], Bayes_estimate[ind]],[min(values), max(values)+.1*(max(values)-min(values))], label = label3[ind], color = 'r')    
    #plt.plot([theta[ind], theta[ind]],[min(values), max(values)+.1*(max(values)-min(values))], label = label2[ind])  
    plt.plot([CI[ind,0], CI[ind,0]],[min(values), max(values)+.1*(max(values)-min(values))], label = r'$95\% \ Credibility \ Interval$', color = 'k', linestyle = '--')
    plt.plot([CI[ind,1], CI[ind,1]],[min(values), max(values)+.1*(max(values)-min(values))], color = 'k', linestyle = '--')    
    plt.ylim([min(values), max(values)+.1*(max(values)-min(values))])
    plt.xlabel(label4[ind])
    plt.ylabel('density')
    plt.xlim([xmin,xmax])
    plt.rc('axes', labelsize=15) 
    plt.legend(loc='best', frameon=False, numpoints=1)
    font = {'size'   : 15}
    plt.rc('font', **font)
    plt.savefig('examples/bastien/Figure/experimental_5_posterior_distribution_'+str(ind)+'.eps', format='eps', dpi=1000) 
np.savetxt('BE.txt', BE)
# Quantify variance for same parameter
lb, ub = [50,5,0.1,0.5e-3,0], [150,20,1.5,3e-3,10]
prior = Uniform(lb = [50,5,0.1,0.5e-3,0],ub = [150,20,1.5,3e-3,10])
pAd, pAg, pT, pF, aT = 98.8, 13.3, 0.495, 0.00129, 7.04
model = Deposition(prior, pAd, pAg, pT, pF, aT, seed = 1)
y_exp = np.load('depo_experimental.npy')
y_obs = model.simulate(60)
np.save('Bayes_estimate_simulation.npy',y_obs)
for ind1 in np.arange(1,5):
    plt.figure()
    plt.plot(y_exp[1:5,0], y_exp[1:5,ind1],'r-*')
    for ind in range(60):
        plt.plot(y_obs[ind][1:5,0], y_obs[ind][1:5,ind1],'k-*')
        
# Generating dataset consisting of 9 patients and mean and std dev
A = np.loadtxt('depo_exp')
T = np.array([0, 20, 60, 120, 300]).reshape(5,1)
B = []
for ind1 in range(11):
    B.append(np.concatenate((T, A[ind1,:].reshape(5,4)), axis=1))
np.save('depo_experimental.npy',B)

# PLot the dataset
label = ['time',r'$\mathcal{S}_{ac}$',r'$\mathbb{N}_{ac}$', r'$\mathbb{N}_{acp}$', r'$\mathbb{N}_{nacp}$']
for index in [1,2,3,4]:        
    plt.figure()
    A = np.zeros(shape=(9,5))
    for ind1 in range(9):
            plt.plot(B[ind1][:,0], B[ind1][:,index], 'k-.')
            A[ind1,:] =  B[ind1][:,index]
    plt.errorbar(B[9][:,0], B[9][:,index], B[10][:,index], marker='^')
    plt.plot(B[9][:,0], np.max(A, axis = 0) , 'r.-')
    plt.plot(B[9][:,0], np.min(A, axis = 0) , 'r.-')
    plt.xlim([0,310])
    plt.xlabel('Time (secs.)')
    plt.ylabel(label[index])   
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #plt.rc('axes', labelsize=15) 
    font = {'size'   : 15}
    plt.rc('font', **font)
    plt.savefig('examples/bastien/Figure/dataset_'+str(index)+'.eps', format='eps', dpi=1000) 