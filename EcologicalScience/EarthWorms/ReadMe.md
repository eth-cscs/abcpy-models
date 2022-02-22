# EarthWorms model (NetLogo)

The models here are implemented with NetLogo (see [here](https://ccl.northwestern.edu/netlogo/)); for now, we wrapped with ABCpy two of the five different models used in the paper [1], in which inference was originally performed in R. Note that these models require an installation of NetLogo and of the Python library `pyNetLogo`. The code was tested with NetLogo 5.3.1 and pyNetLogo 0.4.2. The NetLogo code was written in fact with NetLogo 5, for which `pyNetLogo` supports versions 5.2 and 5.3.    


In the model definition, you need to link to the correct NetLogo path; change line 7 in `model.py`:
  
```netlogo_home = "/path/to/netlogo"```

to your path.

The `GUIDE.pdf` file was provided with the original code and gives more details.

The `src/models` folder contains the original NetLogo model implementation. In `src/R`, the files used for reproducing results in [1] are given for reference. They are however not used by the Python interface, which uses directly the NetLogo ones.


[1] van der Vaart, Elske, et al. "Calibration and evaluation of individual-based models using Approximate Bayesian Computation." Ecological Modelling 312 (2015): 182-190.