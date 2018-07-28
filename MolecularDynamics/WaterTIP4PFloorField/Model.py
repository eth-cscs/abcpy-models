from abcpy.probabilisticmodels import ProbabilisticModel, Continuous, Hyperparameter, InputConnector
import numpy as np
import subprocess, os
import string
from scipy import stats
import time


class TIP4PGromacsOOOH(ProbabilisticModel, Continuous):
    def __init__(self, parameters, name='TIP4PGromacsOOOH'):
        """
        This class implements a probabilistic model simulating rdf_OO, rdf_OH and msd from TIP4P model of water using Gromacs
        with parameters sigma, epsilon and nstep of times.

        Parameters
        ----------
        parameters: list
            Contains the probabilistic models and hyperparameters from which the model derives. Note that the values of the list is not allowed to be smaller than 0.

        name: string
            The name that should be given to the probabilistic model in the journal file.
        """

        if not isinstance(parameters, list):
            raise TypeError('Input of TIP4PGromacsOOOH model is of type list')

        if len(parameters) != 3:
            raise RuntimeError('Input list must be of length 3, containing [sigma, epsilon, nstep].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 3:
            raise ValueError('Number of parameters of TIP4PGromacsOOOH model must be 3.')

        # Check whether input is from correct domain
        sigma = input_values[0]
        epsilon = input_values[1]
        nstep = input_values[2]

        if not isinstance(sigma, (float, np.float64, np.float32, np.float16)) or sigma<=0 :
            return False
        if not isinstance(epsilon, (float, np.float64, np.float32, np.float16)) or epsilon<=0 :
            return False
        if not isinstance(epsilon, (int, np.int64, np.int32, np.int16)):
            return False
        return True

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        # Extract the input parameters
        sigma = input_values[0]
        epsilon = input_values[1]
        nstep = input_values[2]

        # Do the actual forward simulation
        vector_of_k_samples = self.simulate_tpi4p_gromacs(sigma, epsilon, nstep, k)
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result

    def _check_output(self, values):
        return True

    def get_output_dimension(self):
        return 1

    def simulate_tpi4p_gromacs(self, sigma, epsilon, nstep, k):

        resultoutput = []

        for indrep in range(k):
            # The minimum value of nstep can be 50000
            #Creating a folder with random name
            import random
            random = '/Result/'+''.join([random.choice(string.ascii_letters + string.digits) for n in range(32)])

            os.system('mkdir '+random)
            os.system('cp -r '+ '/initial/initial1/.' +' '+random )

            # Write nstep in the correct file
            f = open(random+"/npt.mdp" , 'r+b')
            # get array of lines
            f_content = f.readlines()
            # get middle line
            #middle_line = len(f_content)/2
            # overwrite middle line
            replacement_line_2 = "nsteps = "+str(nstep)+"          ;"+"\n"
            f_content[5] = replacement_line_2.encode('utf-8')
            # return pointer to top of file so we can re-write the content with replaced string
            f.seek(0)
            # clear file content
            f.truncate()
            # re-write the content with the updated content
            f.write(b''.join(f_content))
            # close file
            f.close()



            # Write sigma and epsilon in the correct file
            f = open(random+"/ffnonbonded.itp" , 'r+b')
            # get array of lines
            f_content = f.readlines()
            # get middle line
            #middle_line = len(f_content)/2
            # overwrite middle line
            replacement_line = 'OW_tip4p     8      16.00    0.0000  A   '+str(sigma)+'  '+str(epsilon)+"\n"
            f_content[65] = replacement_line.encode('utf-8')
            # return pointer to top of file so we can re-write the content with replaced string
            f.seek(0)
            # clear file content
            f.truncate()
            # re-write the content with the updated content
            f.write(b''.join(f_content))
            # close file
            f.close()

            # System call to run the gromacs simulation
            stime = time.time()
            p = subprocess.Popen('sh run.sh', shell=True, cwd=random+'/')
            while p.poll() is None and time.time()-stime<3600:
              time.sleep(2)
            os.system('cd ..')

            # Read Output
            rdf_OH = np.matrix(np.loadtxt(random + '/rdf_OH.xvg', skiprows=25))
            rdf_OO = np.matrix(np.loadtxt(random + '/rdf_OO.xvg', skiprows=25))
            msd = np.array(np.loadtxt(random + '/msd.xvg', skiprows=17))

            # Mean of rdf OH
            OH_mean = np.mean(rdf_OH[:, 1])

            # Analyze rdf OH
            rdf_OH = rdf_OH[np.max([np.where(rdf_OH[:, 0] > .1)[0].min(), np.where(rdf_OH[:, 1] == 0)[0].max()]):, :]
            rdf_OH_grad = np.gradient(rdf_OH[:, 1].reshape(-1, ).tolist()[0], rdf_OH[:, 0].reshape(-1, ).tolist()[0])
            asign_OH_grad = np.sign(rdf_OH_grad)
            signchange_OH = np.where(((np.roll(asign_OH_grad, 1) - asign_OH_grad) != 0).astype(int) == 1)[0]

            # Hydrogen bond distance
            OH_min_1 = rdf_OH[signchange_OH[2], 0]
            # Number of hydrogen bonds
            no_H_bonds = np.trapz(rdf_OH[:signchange_OH[2], 1].reshape(-1, ).tolist()[0],
                                  rdf_OH[:signchange_OH[2], 0].reshape(-1, ).tolist()[0])
            # Second minima
            #OH_min_2 = rdf_OH[signchange_OH[4], 0]

            # Analyze rdf OO
            rdf_OO = rdf_OO[np.where(rdf_OO[:, 1] == 0)[0].max():, :]
            rdf_OO_grad = np.gradient(rdf_OO[:, 1].reshape(-1, ).tolist()[0], rdf_OO[:, 0].reshape(-1, ).tolist()[0])
            asign_OO_grad = np.sign(rdf_OO_grad)
            signchange_OO = np.where(((np.roll(asign_OO_grad, 1) - asign_OO_grad) != 0).astype(int) == 1)[0]

            # Mean of rdf OO
            OO_mean = np.mean(rdf_OO[:, 1])

            # First maximum of OO
            OO_max_1 = rdf_OO[signchange_OO[1], 1]
            OO_max_x = rdf_OO[signchange_OO[1], 0]
            # Hydrogen bond distance
            OO_min_1 = rdf_OO[signchange_OO[2], 0]
            # Number of hydrogen bonds
            no_first_nbs = np.trapz(rdf_OO[:signchange_OO[2], 1].reshape(-1, ).tolist()[0],
                                    rdf_OO[:signchange_OO[2], 0].reshape(-1, ).tolist()[0])

            # Analyze MSD
            diff_coeff = stats.linregress(msd[2:, :])[0]

            result = np.array([OH_mean, OH_min_1, no_H_bonds, OO_max_1, OO_max_x, OO_mean, OO_min_1, no_first_nbs, diff_coeff]).reshape(
                1, 9)

            # Delete the random folder
            os.system('rm -r -f ' + random)

            # Save result for output
            resultoutput.append(result)

        return resultoutput

class TIP4PGromacsOO(ProbabilisticModel, Continuous):
    def __init__(self, parameters, name='TIP4PGromacsOO'):
        """
        This class implements a probabilistic model simulating rdf_OO and msd from TIP4P model of water using Gromacs
        with parameters sigma, epsilon and nstep of times.

        Parameters
        ----------
        parameters: list
            Contains the probabilistic models and hyperparameters from which the model derives. Note that the values of the list is not allowed to be smaller than 0.

        name: string
            The name that should be given to the probabilistic model in the journal file.
        """

        if not isinstance(parameters, list):
            raise TypeError('Input of TIP4PGromacsOO model is of type list')

        if len(parameters) != 3:
            raise RuntimeError('Input list must be of length 3, containing [sigma, epsilon, nstep].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 3:
            raise ValueError('Number of parameters of TIP4PGromacsOOOH model must be 3.')

        # Check whether input is from correct domain
        sigma = input_values[0]
        epsilon = input_values[1]
        nstep = input_values[2]

        if not isinstance(sigma, (float, np.float64, np.float32, np.float16)) or sigma<=0 :
            return False
        if not isinstance(epsilon, (float, np.float64, np.float32, np.float16)) or epsilon<=0 :
            return False
        if not isinstance(epsilon, (int, np.int64, np.int32, np.int16)):
            return False
        return True

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        # Extract the input parameters
        sigma = input_values[0]
        epsilon = input_values[1]
        nstep = input_values[2]

        # Do the actual forward simulation
        vector_of_k_samples = self.simulate_tpi4p_gromacs(sigma, epsilon, nstep, k)
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result

    def _check_output(self, values):
        return True

    def get_output_dimension(self):
        return 1

    def simulate_tpi4p_gromacs(self, sigma, epsilon, nstep, k):

        resultoutput = []

        for indrep in range(k):
            # The minimum value of nstep can be 50000
            #Creating a folder with random name
            import random
            random = '/Result/'+''.join([random.choice(string.ascii_letters + string.digits) for n in range(32)])

            os.system('mkdir '+random)
            os.system('cp -r '+ '/initial/initial1/.' +' '+random )

            # Write nstep in the correct file
            f = open(random+"/npt.mdp" , 'r+b')
            # get array of lines
            f_content = f.readlines()
            # get middle line
            #middle_line = len(f_content)/2
            # overwrite middle line
            replacement_line_2 = "nsteps = "+str(nstep)+"          ;"+"\n"
            f_content[5] = replacement_line_2.encode('utf-8')
            # return pointer to top of file so we can re-write the content with replaced string
            f.seek(0)
            # clear file content
            f.truncate()
            # re-write the content with the updated content
            f.write(b''.join(f_content))
            # close file
            f.close()



            # Write sigma and epsilon in the correct file
            f = open(random+"/ffnonbonded.itp" , 'r+b')
            # get array of lines
            f_content = f.readlines()
            # get middle line
            #middle_line = len(f_content)/2
            # overwrite middle line
            replacement_line = 'OW_tip4p     8      16.00    0.0000  A   '+str(sigma)+'  '+str(epsilon)+"\n"
            f_content[65] = replacement_line.encode('utf-8')
            # return pointer to top of file so we can re-write the content with replaced string
            f.seek(0)
            # clear file content
            f.truncate()
            # re-write the content with the updated content
            f.write(b''.join(f_content))
            # close file
            f.close()

            # System call to run the gromacs simulation
            p = subprocess.Popen('sh run.sh', shell=True, cwd=random+'/')
            p.wait()
            os.system('cd ..')

            # Read Output
            rdf_OH = np.matrix(np.loadtxt(random + '/rdf_OH.xvg', skiprows=25))
            rdf_OO = np.matrix(np.loadtxt(random + '/rdf_OO.xvg', skiprows=25))
            msd = np.array(np.loadtxt(random + '/msd.xvg', skiprows=17))

            # Mean of rdf OH
            OH_mean = np.mean(rdf_OH[:, 1])

            # Analyze rdf OH
            rdf_OH = rdf_OH[np.max([np.where(rdf_OH[:, 0] > .1)[0].min(), np.where(rdf_OH[:, 1] == 0)[0].max()]):, :]
            rdf_OH_grad = np.gradient(rdf_OH[:, 1].reshape(-1, ).tolist()[0], rdf_OH[:, 0].reshape(-1, ).tolist()[0])
            asign_OH_grad = np.sign(rdf_OH_grad)
            signchange_OH = np.where(((np.roll(asign_OH_grad, 1) - asign_OH_grad) != 0).astype(int) == 1)[0]

            # Hydrogen bond distance
            OH_min_1 = rdf_OH[signchange_OH[2], 0]
            # Number of hydrogen bonds
            no_H_bonds = np.trapz(rdf_OH[:signchange_OH[2], 1].reshape(-1, ).tolist()[0],
                                  rdf_OH[:signchange_OH[2], 0].reshape(-1, ).tolist()[0])
            # Second minima
            #OH_min_2 = rdf_OH[signchange_OH[4], 0]

            # Analyze rdf OO
            rdf_OO = rdf_OO[np.where(rdf_OO[:, 1] == 0)[0].max():, :]
            rdf_OO_grad = np.gradient(rdf_OO[:, 1].reshape(-1, ).tolist()[0], rdf_OO[:, 0].reshape(-1, ).tolist()[0])
            asign_OO_grad = np.sign(rdf_OO_grad)
            signchange_OO = np.where(((np.roll(asign_OO_grad, 1) - asign_OO_grad) != 0).astype(int) == 1)[0]

            # Mean of rdf OO
            OO_mean = np.mean(rdf_OO[:, 1])

            # First maximum of OO
            OO_max_1 = rdf_OO[signchange_OO[1], 1]
            OO_max_x = rdf_OO[signchange_OO[1], 0]
            # Hydrogen bond distance
            OO_min_1 = rdf_OO[signchange_OO[2], 0]
            # Number of hydrogen bonds
            no_first_nbs = np.trapz(rdf_OO[:signchange_OO[2], 1].reshape(-1, ).tolist()[0],
                                    rdf_OO[:signchange_OO[2], 0].reshape(-1, ).tolist()[0])

            # Analyze MSD
            diff_coeff = stats.linregress(msd[2:, :])[0]

            result = np.array([OO_max_1, OO_max_x, OO_mean, OO_min_1, no_first_nbs, diff_coeff]).reshape(
                1, 6)

            # Delete the random folder
            os.system('rm -r -f ' + random)

            # Save result for output
            resultoutput.append(result)

        return resultoutput


