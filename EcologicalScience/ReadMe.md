# Ecological Science

We give here some models from ecological science.

- [Ricker model](https://en.wikipedia.org/wiki/Ricker_model) is a discrete population model; we use the version described in [1] which adds some stochastic noise
- [Lotka-Volterra model](https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations) is the standard prey-predator dynamic model, defined in terms of a coupled ODE; we implement the version described in Appendix T.10 in [2], where the ODE is integrated and observations are made at 10 evenly-spaced time points inside the integration interval. Both the number of preys and predators are observed, with an (optional) lognormal noise.


[1] S. N. Wood. Statistical inference for noisy nonlinear ecological dynamic systems. Nature, 466(7310):1102–1104, Aug. 2010.

[2] Lueckmann, Jan-Matthis, et al. "Benchmarking Simulation-Based Inference." arXiv preprint arXiv:2101.04653 (2021), url https://arxiv.org/abs/2101.04653.