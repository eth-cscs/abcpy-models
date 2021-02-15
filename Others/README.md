# Other models

We have here some other models: 


- Iid models (see details in the folder)
- Timeseries models (see details in the folder)
- Bivariate normal parametrized by `[mu1, mu2, sigma1, sigma2, rho]`
- Multivariate g-and-k (see more info below)
- Mixture of two bivariate gaussian components, with fixed covariance matrices and parametrized by `[p, mu0_0, mu0_1, mu1_0, mu1_1]`, `p` being the mixture weight of the component with mean `mu1`.


## Multivariate g-and-k

The definition of this follows Section 4.5 in [1]; the code implementation was inspired by [2]. Drawing from this distribution is cheap, while evaluating the likelihood is costly as you need to numerically invert an equation. We provide code for doing it anyway. 

Please find more details on multivariate quantile distributions and the computation of the likelihood in [3].


The `g-and-k.py` code was inspired by the R package `gk`  (described in [2], source code https://github.com/dennisprangle/gk) which was made available under the GNU General Public License v2.0; the code therein is therefore made available under the same license, which is reported for ease of use in `LICENSE_g_and_k.md`

[1] Jiang, Bai. "Approximate Bayesian computation with Kullback-Leibler divergence as data discrepancy."
International Conference on Artificial Intelligence and Statistics. PMLR, 2018.

[2] Prangle, Dennis. "gk: An R Package for the g-and-k and Generalised g-and-h Distributions."
arXiv preprint arXiv:1706.06889 (2017), url https://arxiv.org/abs/1706.06889.

[3] Drovandi, Christopher C., and Anthony N. Pettitt. "Likelihood-free Bayesian estimation of multivariate quantile distributions." Computational Statistics & Data Analysis 55.9 (2011): 2541-2556.