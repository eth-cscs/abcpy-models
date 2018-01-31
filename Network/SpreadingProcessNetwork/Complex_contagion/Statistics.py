import numpy as np
from abcpy.statistics import Statistics
import networkx as nx

class Diffusion_Identity_statistics(Statistics):
    """
    This class implements statistics for linear diffusion.
    """
    def __init__(self, degree = 1, cross = False):
        self.degree = degree
        self.cross = cross
        

    def statistics(self, data):  
        
       return data           
        

