# ARMA model

We give here an ARMA model which can be set up with different lags for both the AR and MA component.
See [here](https://en.wikipedia.org/wiki/Autoregressive%E2%80%93moving-average_model#ARMA_model) for more information.
Under the hood, it uses the `statsmodels` package to implement that.

As possible set of statistics, we give the Autocorrelation. You can choose up to which order to use.

Notice that the exact likelihood of ARMA models can be computed in the case in which only the AR or MA components are used independently. We give in `likelihoods.py` the likelihood for AR(2) and MA(2) case.  