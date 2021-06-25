# Bass model (NetLogo)

The model here are implemented with NetLogo (see [here](https://ccl.northwestern.edu/netlogo/)); for now, we wrapped with ABCpy the model used in the paper [1]. Note that these models require an installation of NetLogo and of the Python library `pyNetLogo`. In the model definition, you need to link to the correct netlogo path; change line 170 in `model.py`:
  
```netlogo_home = "/path/to/netlogo"```

to your path.



[1] Watson, Joseph W, et al. "Incorporating environmental variability in a spatially-explicit individual-based model of European sea bass".