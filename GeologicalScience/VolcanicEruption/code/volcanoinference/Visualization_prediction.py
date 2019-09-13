import numpy as np
import pylab as plt
from distance_learning.utilities import read_simulations

a, b, output = read_simulations('output1.txt')
a, b, output = read_simulations('predicted_data.txt')

print(output.shape)

obs_data = [
    np.array([70.0, 100.0, 170.0, 150.0, 290.0, 220.0, 235.0, 180.0, 275.0, 175.0, 170.0,
              165.0, 140.0, 100.0, 100.0, 120.0, 110.0, 80.0, 25.0, 20.0, 30.0, 40.0, 30.0,
              20.0, 15.0, 30.0, 70.0, 70.0, 100.0, 70.0, 75.0, 75.0, 100.0, 110.0, 115.0, 80.0,
              125.0, 280.0, 350.0, 350.0, 180.0, 25.0, 40.0, 40.0, 55.0, 70.0, 175.0, 260.0,
              200.0, 140.0, 50.0, 10.0, 120.0, 300.0, 250.0, 60.0, 90.0, 105.0, 150.0, 170.0,
              120.0, 110.0, 120.0, 190.0, 190.0, 130.0, 130.0, 110.0, 70.0, 90.0, 85.0, 65.0]).reshape(-1, )]

mean_prediction = np.mean(output, axis=0)
quantile_prediction1 = np.quantile(output, .975, axis=0)
quantile_prediction2 = np.quantile(output, .025, axis=0)
min_prediction = np.min(output, axis=0)
max_prediction = np.max(output, axis=0)

# plt.figure()
fig, ax = plt.subplots()
ppdata = ax.plot(obs_data[0], 'k*')
ppmean = ax.plot(mean_prediction, 'r*')
bp = ax.boxplot(output, notch=True, meanline=True, showfliers=False, patch_artist=True,
                boxprops=dict(facecolor=(0, 0, 0, 0)))
ax.legend([ppdata[0], bp["boxes"][0], ppmean[0]],
          ['2450BP Pululagua deposit', 'Prediction uncertainty', 'Mean Prediction'])
# ax = plt.fill_between(range(72), quantile_prediction1, quantile_prediction2, alpha=0.2, label='95% C.I. of Predictive Distribution')
# plt.plot(quantile_prediction1, 'b-')
# plt.plot(quantile_prediction2, 'b-')
ax.set_xticks([0, 20, 40, 60])
ax.set_xticklabels([0, 20, 40, 60])
# ax.legend()
ax.set_ylabel('Tephra deposit ' + r'$[kg/m^2]$')
ax.set_xlabel('Location')
fig.savefig('figures/PredictionCheck.eps')
plt.show()
