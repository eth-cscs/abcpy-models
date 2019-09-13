from abcpy.probabilisticmodels import ProbabilisticModel, InputConnector
from model import model
import os
import uuid
from fileio import *
import numpy as np

class Volcano(ProbabilisticModel):

    def __init__(self, parameters, name='Volcano'):
        # We expect input of type parameters = [U0, L]
        if not isinstance(parameters, list):
            raise TypeError('Input of volcano model is of type list')

        if len(parameters) != 2:
            raise RuntimeError('Input list must be of length 2, containing [U0, L]')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)


    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 2:
            raise ValueError('Number of parameters is 2')
        return True


    def _check_output(self, values):
        if not isinstance(values, np.ndarray):
            raise ValueError('Output of the model is always should be a numpy array')

        if values.shape[0] != 72:
            raise ValueError('Output shape should be of dimension 72')

        return True

    def get_output_dimension(self):
        return 72

    def forward_simulate(self, input_values, k, rng=np.random.RandomState(), mpi_comm=None):
        # do i k times
        results = []
        for i in range(k):
            # determine a file name based on parameters and a unique id
            if mpi_comm.Get_rank() == 0:
                id = uuid.uuid4()
            else:
                id = None
            id = mpi_comm.bcast(id)
            U0 = input_values[0]
            L = input_values[1]
            seed = 341 #rng.randint(1, 10000)
            filename = "deposit-U0-"+str(U0)+"-L-"+str(L)+"-"+str(id)+".h5"

            # run the model
            success = model(mpi_comm, U0, L, seed, filename)

            # extract statistics from generated data
            if mpi_comm.Get_rank() == 0:
                if success == 1:
                    results.append(self.extract_deposition_vector(filename))
                    os.remove(filename)
                else:
                    results.append(np.random.uniform(0, 1, 72))
                    
        #print('Results: '+str(results))
        result = None
        if mpi_comm.Get_rank()==0:
            result = [np.array([results[i]]).reshape(-1,) for i in range(k)]
        #print('Result: '+str(result))
        f = open("output.txt", "a")
        f.write(str(U0)+" "+str(L)+" "+str(result)+"\n")
        f.close()
        result = mpi_comm.bcast(result)
        return result

    # extract the deposition vector from given deposition file
    def extract_deposition_vector(self, filename):
        measures = readPhysicalParam("Physical Parameters BF2.csv")
        points = massareaPoints(filename, measures, "all", 3.0)
        measures = array(list(zip(*points))[0])
        compute = array(list(zip(*points))[1])
        return compute
