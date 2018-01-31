import numpy as np
from abcpy.models import Model
from examples.bastien.param_5.deposition_model import deposition_model


class Deposition(Model):
    def __init__(self, prior, pAd = None, pAg = None, pT = None, pF = None, aT = None, seed=None):
        """
        Parameters
        ----------
        prior: abcpy.distributions.Distribution
            Prior distribution
        pAd: integer
            pAdhesion
        pAg: float
            pAggreg
        pT: float
            pTop
        pF: float
            pfilling; adhesion of albumine
        seed: int, optional
            Initial seed. The default value is generated randomly.             
        """    
        self.prior = prior
        # test provided model parameters
        if (pAd == None) + (pAg == None) + (pT == None) + (pF == None) + (aT == None)<4 and (pAd == None) + (pAg == None) + (pT == None) + (pF == None) + (aT == None)>0:
            raise ValueError("All or neither of the model parameters have to be provided.")

        # set model parameters directly if specified
        if pAd != None and pAg != None and pT != None and pF != None and aT != None:
            if self.set_parameters(np.array([pAd,pAg,pT,pF, aT])) == False:
                raise ValueError("The parameter values are out of the model parameter domain.")
        else:
            self.sample_from_prior()
        self.rng = np.random.RandomState(seed)
        

    def set_parameters(self, theta):
        theta = np.array(theta)
        if theta.shape[0] > 5: return False
        if theta[0] <= 50 or theta[0] >= 150: return False
        if theta[1] <= 5 or theta[1] >= 20: return False
        if theta[2] <= 0.1 or theta[2] >= 1.5: return False
        if theta[3] <= 0.5e-3  or theta[3] >= 3e-3: return False        
        if theta[4] <= 0  or theta[4] >= 10: return False 

        self.pAd = theta[0]
        self.pAg = theta[1]
        self.pT = theta[2]
        self.pF = theta[3]
        self.aT = theta[4]
        return True

    def get_parameters(self):
        return np.array([self.pAd, self.pAg, self.pT, self.pF, self.aT])

    def sample_from_prior(self):
        sample = self.prior.sample(1).reshape(-1)
        self.set_parameters(sample)

    def simulate(self, k):
        seed = self.rng.randint(np.iinfo(np.int32).max)

        nrow = 5
        ncol = 5
        mshape = ncol*nrow
        rshape = nrow*ncol*k
        results = np.reshape(deposition_model(rshape, k, self.pAd, self.pAg, self.pT, self.pF, self.aT, seed),[k,mshape])
        result = [None]*k
        for ind in range(k):
            result[ind] = np.reshape(results[ind],[nrow,ncol])
	
        return result
