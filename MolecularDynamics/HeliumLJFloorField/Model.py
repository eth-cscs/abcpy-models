from abcpy.probabilisticmodels import ProbabilisticModel, Continuous, InputConnector
import numpy as np
import subprocess, os, errno
import string
import time
from pathlib import Path

class Helium(ProbabilisticModel, Continuous):
    """
    This class is an re-implementation of the `abcpy.continousmodels.Normal` for documentation purposes.
    """

    def __init__(self, parameters, name='Helium'):
        # We expect input of type parameters = [mu, sigma]
        if not isinstance(parameters, list):
            raise TypeError('Input of Helium model is of type list')

        if len(parameters) != 2:
            raise RuntimeError('Input list must be of length 2, containing [sigma, epsilon].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)


    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 2:
            raise ValueError('Number of parameters of IterMap model must be 2.')

        # Check whether input is from correct domain
        sigma = input_values[0]
        epsilon = input_values[1]

        if sigma <=0 or sigma >= 1 or epsilon <=0 or epsilon >= 1:
            return False

        return True


    def _check_output(self, values):
        return True


    def get_output_dimension(self):
        return 1


    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        # Extract the input parameters
        sigma = input_values[0]
        epsilon = input_values[1]

        # Do the actual forward simulation
        vector_of_k_samples = self.heliumlammps(sigma, epsilon, k)
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result

    def heliumlammps(self, sigma, epsilon, k):
        """
        Description of the Function

        :param sigma: float (should be )
        Description of epsilon
        :param epsilon: float (should be greater than zero)
        Description of epsilon
        :param k: int
        Number of data to be simulated
        :return: list
        This returns a list with k elements, where each element is a numpy array consisting of a time-series
        """
        ###############################################################################################
        ###############################################################################################
        # Creating a folder with random name
        import random
        import sys

        result = []
        for ind_k in range(k):
            random = '/Result/'+''.join([random.choice(string.ascii_letters + string.digits) for n in range(32)])
            os.system('mkdir ' + random)
            os.system('cp -r ' + '/Helium/.' + ' ' + random)

            # Write sigma and epsilon in the correct file
            f = open(random + "/helium.in", 'r+b')
            # get array of lines
            f_content = f.readlines()
            # get middle line
            # overwrite middle line
            replacement_line_1 = "variable                  sigma equal " + str(sigma) + " " + "# nm" + "\n"
            replacement_line_2 = "variable                  eps equal " + str(epsilon) + " " + "# (ag*nm^2)/(ns^2*K)" + "\n"
            f_content[27] = replacement_line_1.encode('utf-8')
            f_content[28] = replacement_line_2.encode('utf-8')
            # return pointer to top of file so we can re-write the content with replaced string
            f.seek(0)
            # clear file content
            f.truncate()
            # re-write the content with the updated content
            f.write(b''.join(f_content))
            # close file
            f.close()

            iterationrun = 0
            run = True
            while run and iterationrun<100:
                try:
                    stime = time.time()
                    ## The user need to fill the address of executable of LAMPS in the following line
                    p = subprocess.run('LAMPS_exec -in helium.in', shell=True,cwd=random + '/', timeout=2000)
                    print(random+' : Program ran successfully. Time taken :'+str(time.time()-stime) + '\n')  
                    if Path(random + '/' + 'log.lammps').is_file() is False:
                        raise ValueError('Output not created')
                    else:
                        run = False
                except OSError as e:
                    print(random+' : Error occurred: '  + str(e.errno) +'. Will try again.' +'\n')
                except subprocess.TimeoutExpired:
                    print(random+' : Process ran too long. Will try again.'+'\n')
                except (ValueError, IndexError):
                    print(random+' : Output not created'+'\n')   
                iterationrun += 1

            # Read Output
            if Path(random + '/' + 'log.lammps').is_file():
                f = open(random + '/' + 'log.lammps', "r")
                lines = f.readlines()[371:1371]
                relentropy = []
                for x in lines:
                    relentropy.append(float(x.split()[6]))
                f.close()
                # successfullyran - make True to break the if loop
                notsuccessfullyran = False
                # Create array of simulated data
                result.append(np.array(relentropy))
            else:
                result.append(np.random.uniform(0,1,1000))

            # Delete the random folder
            os.system('rm -r -f ' + random)
        return result
