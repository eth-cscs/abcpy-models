# Lorenz 95 Weather prediction model
The Lorenz model is a modification of the original weather prediction model of Lorenz (1995)
when fast climate variables are unobserved (Wilks (2005)). To perform inference using ABC, we provide
summary statistics suggested by Hakkarainen et al. (2012).

## Reference
1. Lorenz E (1995). “Predictability: a problem partly solved.” In Proceedings of the Seminar on
Predictability, 4-8 September 1995, volume 1, pp. 1–18. European Center on Medium Range
Weather Forecasting, European Center on Medium Range Weather Forecasting, Shinfield
Park, Reading.
2. Wilks DS (2005). “Effects of stochastic parametrizations in the Lorenz ’96 system.” Quarterly
Journal of the Royal Meteorological Society, 131(606), 389–407. doi:10.1256/qj.04.03.
3. Hakkarainen J, Ilin A, Solonen A, Laine M, Haario H, Tamminen J, Oja E, Järvinen H (2012).
“On closure parameter estimation in chaotic systems.” Nonlinear Processes in Geophysics,
19(1), 127–143. doi:10.5194/npg-19-127-2012.

## Dependencies 
- abcpy0.5.1

## Setup
You can run an example ABC inference scheme on Lorenz model by 'python3 LorenzInference.py'

