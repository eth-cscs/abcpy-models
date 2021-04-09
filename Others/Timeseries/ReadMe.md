# ARMA model

We give here an ARMA model which can be set up with different lags for both the AR and MA component.
See [here](https://en.wikipedia.org/wiki/Autoregressive%E2%80%93moving-average_model#ARMA_model) for more information.
Under the hood, it uses the `statsmodels` package to implement that.

As possible set of statistics, we give the Autocorrelation. You can choose up to which order to use.

Notice that the exact likelihood of ARMA models can be computed in the case in which only the AR or MA components are used independently. We give in `likelihoods.py` the likelihood for AR(2) and MA(2) case.  

This implementation requires `statsmodels==0.9.0`; you can install it with: 

    pip install statsmodels==0.9.0


In some works, the MA2 model is used with a uniform prior on a triangular region for [theta1, theta2]. We provide a reparametrized model
which takes as inputs parameters [R1,R2] which, if they are given a Uniform[0,1] prior, result in theta1, theta2 being distributed uniformly on the correct 
triangular region (see for instance in [1])

[1] Marin, J. M., Pudlo, P., Robert, C. P., & Ryder, R. J. (2012). Approximate Bayesian computational methods. 
Statistics and Computing, 22(6), 1167-1180.