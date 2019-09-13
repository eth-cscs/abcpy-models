from abcpy.output import Journal
import numpy as np
import pylab as plt
from scipy.stats import gaussian_kde

type = 'simulated'
# type = 'observed'

if type == 'simulated':
    #filename = 'VolcanojournalAPMCABC_triplet_simulated'
    filename = 'VolcanojournalAPMCABC_simulated'
else:
    #filename = 'VolcanojournalAPMCABC_triplet-fromfile_Pululagua'
    filename = 'VolcanojournalAPMCABC_pululagua'

journal = Journal.fromFile(filename + '.jrnl')

# for ind in range(len(journal.weights)):
#     journal.weights[ind] = journal.weights[ind]/sum(journal.weights[ind])
# journal.save(filename + '_normalized.jrnl')

mean = journal.posterior_mean() 
print(journal.posterior_mean())
print(journal.posterior_cov())
k = -1
prioru = np.concatenate(journal.get_parameters(0)['U0'])
priorl = np.concatenate(journal.get_parameters(0)['L'])
postu = np.concatenate(journal.get_parameters(k)['U0']).reshape(-1,)
postl = np.concatenate(journal.get_parameters(k)['L']).reshape(-1,)
weights = np.concatenate(journal.get_weights(k))


#weight = np.concatenate(journal.get_weights(k)).reshape(-1,)
#Indices = np.random.choice(postu.shape[0], size=500, replace=True, p = weight/sum(weight))
#postu = postu[Indices]
#postl = postl[Indices]

# Peform the kernel density estimate
xmin, xmax = 100, 300
ymin, ymax = 30, 100
#xmin, xmax = np.min(postu), np.max(postu)
#ymin, ymax = np.min(postl), np.max(postl)
X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([postu, postl])
kernel = gaussian_kde(values, weights=weights, bw_method=.8)
Z = np.reshape(kernel(positions).T, X.shape)
Z = Z/np.max(Z)


fig = plt.figure()
ax = fig.gca()
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
#ax.plot(postu, postl, 'r*', markersize=.5, label='Posterior Samples')
ax.plot(mean['U0'], mean['L'], 'k+', markersize=20, label='Posterior Mean')
#ax.plot(np.mean(postu), np.mean(postl), 'k+', markersize=20, label='Posterior Mean')
# PLot true value for simulated data
if type == 'simulated':
    plt.plot(173.87, 84.55, 'b+', markersize=20, label='True Value')
# Contourf plot
cfset = ax.contourf(X, Y, Z, cmap='Blues')
## Or kernel density estimate plot instead of the contourf plot
#ax.imshow(np.rot90(f), cmap='Blues', extent=[xmin, xmax, ymin, ymax])
# Contour plot
cset = ax.contour(X, Y, Z, levels=[.1, .5, .9], colors='k')
# Label plot
ax.clabel(cset, inline=1, fontsize=10)
ax.set_xlabel(r'$U_0\ [m/s]$', fontsize=15)
ax.set_ylabel(r'$R_0\ [m]$', fontsize=15)
ax.legend(loc='best', frameon=False, numpoints=1, fontsize=15)
plt.tight_layout()
plt.savefig('figures/' + filename + '-no_post_samples.eps', format='eps', dpi=1000)
plt.show()



