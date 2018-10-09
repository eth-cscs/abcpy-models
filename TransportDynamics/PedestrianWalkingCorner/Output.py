import numpy as np
from abcpy.output import Journal
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# Reading your journal file
claudio_journal = Journal.fromFile('apmcabc_fakeobs1.jrnl')

# Read what were the parameters used for this specific inference scheme which produced this result
# each key : value of the key
print('Lets find out the configuration of the inference: '+str(claudio_journal.configuration)+'\n')
# It shows you have run APMCABC for 4 steps with 10000 n_samples and etc.
# Most important is the epsilon_arr - This tells me the 'threshold vaklue' used at each of the 4 steps, this is automatically chosen.
# In this case, your final epsilon value is 0.1696, which is rather large.
# Our goal is to take this value as close to zeor possible
# TIPS: You have to choose big steps more than 4 and see how many do you need to lets say bring them down to 0.01?
# One thing you can do: start the inference APMCABC using the final samples stored in your previous journal file, which is
# in this case 'apmcabc_fakeobs1.jrnl'. To that in the inputs of APMCABC sample put journal_file = 'apmcabc_fakeobs1.jrnl'


# The posterior distribution is approximated by the 10000 (as n_samples=10000) samples drawn from the approximate posterior distribution
# First we read the posterior samples in samples_dictionary which is a python dictionary object
# (Here the samples we consider are stored at the end meaning after step 4)
samples_dictionary = claudio_journal.get_parameters()
# We first check what were parameters you were infering, this is stored as the keys
print('The parameters for which inference was done: '+str(samples_dictionary.keys())+'\n')
# We also get the weights corresponding to each sample
weights = np.array(claudio_journal.get_weights()).reshape(10000,)
# Normalize weights
weights = weights/sum(weights)
print('The corresponding weights: '+str(weights.shape)+'\n')

# Samples corresponding to the parameters 'bufferAngle' are extracted, converted to an numpy array and then resampled
# using the weight
samples_bufferAngle = np.array(samples_dictionary['bufferAngle'])
# The following step is resampling step, which is crucuial
samples_bufferAngle = samples_bufferAngle[np.random.choice(10000, 10000, p=weights),0]
print('How many samples?: '+str(samples_bufferAngle.shape)+'\n')
# Now we plot marginal posterior distribution of bufferAngle using the above samples
# We just use a histrogram plot here
# the histogram of the data
plt.figure()
n, bins, patches = plt.hist(samples_bufferAngle, 50, density=True, facecolor='g', alpha=0.75)
plt.xlabel('bufferAngle')
plt.ylabel('Probability')
plt.title('Marginal posterior distribution of bufferAngle')
# Saved as an eps file, which you can check
plt.savefig('bufferAngle_marginal.eps')

# Lets do the same for kD
samples_kD = np.array(samples_dictionary['kD'])
samples_kD = samples_kD[np.random.choice(10000, 10000, p=weights),0]
plt.figure()
n, bins, patches = plt.hist(samples_kD, 50, density=True, facecolor='g', alpha=0.75)
plt.xlabel('kD')
plt.ylabel('Probability')
plt.title('Marginal posterior distribution of kD')
plt.savefig('kD_marginal.eps')

# Lets do the same for decay
samples_decay = np.array(samples_dictionary['decay'])
samples_decay = samples_decay[np.random.choice(10000, 10000, p=weights),0]
print(samples_decay.shape)
plt.figure()
n, bins, patches = plt.hist(samples_decay, 50, density=True, facecolor='g', alpha=0.75)
plt.xlabel('decay')
plt.ylabel('Probability')
plt.title('Marginal posterior distribution of decay')
plt.savefig('decay_marginal.eps')


# Now lets plot the joint marginal posterior distribution of kD and decay (You can similarly plot for any other pair)
# Here I will use a Gaussian smoothing kernel, which we could have used above in place of histogram
xmin, xmax = min(samples_kD), max(samples_kD)
ymin, ymax = min(samples_decay), max(samples_decay)
X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([samples_kD, samples_decay])
kernel = gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)
plt.figure()
plt.plot(samples_kD,samples_decay,'.',markersize=1,color = '0.8',alpha = 0.2)
CS = plt.contour(X, Y, Z, 6, colors = 'k',linestyles='solid')
plt.clabel(CS, fontsize=9, inline=1)
plt.xlabel(r'$kD$', fontsize=30)
plt.ylabel(r'$decay$', fontsize=30)
plt.tight_layout()
plt.legend(loc='upper left', frameon=False, numpoints=1,fontsize=25)
plt.savefig('kD_decay_marginal.eps', format='eps', dpi=1000)