# EarthWorms model (NetLogo)

The models here are implemented with NetLogo (see [here](https://ccl.northwestern.edu/netlogo/)); for now, we wrapped with ABCpy two of the five different models used in the paper [1], which inference was originally performed in R. Note that these models require an installation of NetLogo and of the Python library `pyNetLogo`. In the model definition, you need to link to the correct netlogo path; change line 7 in `model.py`:
  
```netlogo_home = "/path/to/netlogo"```

to your path.

The `GUIDE.pdf` file was provided with the original code and gives more details.



[1] van der Vaart, Elske, et al. "Calibration and evaluation of individual-based models using Approximate Bayesian Computation." Ecological Modelling 312 (2015): 182-190.