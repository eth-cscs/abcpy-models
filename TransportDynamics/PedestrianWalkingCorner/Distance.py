from abcpy.distances import Distance
import copy
import math
import numpy as np


class Absolute(Distance):
    """
    This distance is simply the Euclidean distance between the heatmaps
    The Euclidean distance between each observation and each simulation is computed
    The final outcome is the average difference for each cell among simulation and observation

    The maximum value of the distance is np.inf.
    """

    def __init__(self):

        # Since the observations do always stay the same, we can save the
        #  summary statistics of them and not recalculate it each time
        self.s1 = None
        self.data_set = None

    def distance(self, d1, d2):
        """Calculates the distance between two datasets.

        Parameters
        ----------
        d1, d2: list
            A list, containing a list describing the data set
        """
        if not isinstance(d1, list):
            raise TypeError('Data is not of allowed types')
        if not isinstance(d2, list):
            raise TypeError('Data is not of allowed types')

        # Observation
        resObservationHeatmap = []
        resObsPos = []
        for repi in range(len(d1)):
            # Extract fundamental information
            nCell = d1[repi][1].astype(int)
            nTime = d1[repi][0].astype(int)
            resObservationHeatmap.append(
                copy.deepcopy(np.array(d1[repi][2:2 + int(pow(nCell, 2))]).reshape(nCell, nCell)))
            resObsPos.append(np.array(d1[repi][2 + int(pow(nCell, 2)) + nTime:]))

        # Simulation
        resSimulationHeatmap = []
        resSimPos = []
        for repi in range(len(d2)):
            # Extract fundamental information
            nCell = d2[repi][1].astype(int)
            nTime = d2[repi][0].astype(int)
            resSimulationHeatmap.append(
                copy.deepcopy(np.array(d2[repi][2:2 + int(pow(nCell, 2))]).reshape(nCell, nCell)))
            resSimPos.append(np.array(d2[repi][2 + int(pow(nCell, 2)) + nTime:]))

        ## ACTUALLY COMPUTE DISTANCE
        # Compute distance using heatmap
        tempDistance = 0
        for simNumber in range(len(d2)):
            for obsNumber in range(simNumber, len(d1)):
                tempDistance += sum(sum(abs(resObservationHeatmap[obsNumber] - resSimulationHeatmap[simNumber]))) / pow(
                    nCell, 2)
                # tempDistance += self.all_dist(resObsPos[obsNumber], resSimPos[simNumber])

        return tempDistance / (len(d1) * len(d2))

    def dist_max(self):
        return np.inf

    def all_dist(self, resObsPos, resSimPos):
        if len(resSimPos) > len(resObsPos):
            resObsPos = resObsPos.tolist() + [0 for i in range(len(resSimPos) - len(resObsPos))]
        else:
            resSimPos = resSimPos.tolist() + [0 for i in range(len(resObsPos) - len(resSimPos))]
        return sum(abs(resObsPos - resSimPos)) / len(resSimPos)

class DistanceType1(Distance):
    """
    This distance is simply the Euclidean distance between the heatmaps
    The Euclidean distance between each observation and each simulation is computed
    The final outcome is the average difference for each cell among simulation and observation

    The maximum value of the distance is np.inf.
    """

    def __init__(self, statistics):
        self.statistics_calc = statistics

        # Since the observations do always stay the same, we can save the
        #  summary statistics of them and not recalculate it each time
        self.s1 = None
        self.data_set = None
        
    def distance(self, d1, d2):
        """Calculates the distance between two datasets.

        Parameters
        ----------
        d1, d2: list
            A list, containing a list describing the data set
        """
        if not isinstance(d1, list):
            raise TypeError('Data is not of allowed types')
        if not isinstance(d2, list):
            raise TypeError('Data is not of allowed types')
        
        # Distance main algorithm
        def single_distance(obsHeatmap,simHeatmap):
            return abs(sum(sum(pow(obsHeatmap - simHeatmap,2))))

        # Extract summary statistics from the dataset
        if (self.s1 is None or (np.asarray(self.data_set) == np.asarray(d1)).all() == False):
            self.s1 = self.statistics_calc.statistics(d1)
            self.data_set = d1
        s2 = self.statistics_calc.statistics(d2)

        ## EXTRACT MATRICES FROM LIST
        obsRep, simRep = len(self.s1), len(s2)

        # Observation
        resObservationHeatmap = []
        for repi in range(obsRep):
            # Extract fundamental information
            nCells = int(self.s1[repi,1])
            # Extract observation heat map
            repObservationHeatmap = np.zeros((nCells,nCells))
            for i in range(nCells):
                index = 2 + i*nCells
                repObservationHeatmap[i,:] = self.s1[repi,index:index+nCells]
            # Add to global quantities
            resObservationHeatmap.append(copy.deepcopy(repObservationHeatmap))
        
        # Simulation
        resSimulationHeatmap = []
        for repi in range(simRep):
            print(s2[repi, :4])
            # Extract fundamental information
            nCells = int(s2[repi,1])
            # Extract simulation heat map
            repSimulationHeatmap = np.zeros((nCells,nCells))
            for i in range(nCells):
                index = 2 + i*nCells
                repSimulationHeatmap[i,:] = s2[repi,index:index+nCells]    
            # Add to global quantities
            resSimulationHeatmap.append(copy.deepcopy(repSimulationHeatmap))
        
        ## ACTUALLY COMPUTE DISTANCE
        # Compute distance using heatmap
        tempDistance = 0
        for simNumber in range(simRep):
            for obsNumber in range(obsRep):
                tempDistance += single_distance(resObservationHeatmap[obsNumber],resSimulationHeatmap[simNumber])
        return tempDistance/(nCells**2 * simRep * obsRep)

    def dist_max(self):
        return np.inf
        
class DistanceType2(Distance):
    """
    This distance is Euclidean distance between the sorted heatmaps
    The idea is that is not so important where the maximum value is but how big it is and how are the number distributed
    This is simply a guess and this distance was defined to check if such idea work or not

    The maximum value of the distance is np.inf.
    """

    def __init__(self, statistics):
        self.statistics_calc = statistics

        # Since the observations do always stay the same, we can save the
        #  summary statistics of them and not recalculate it each time
        self.s1 = None
        self.data_set = None
        
    def distance(self, d1, d2):
        """Calculates the distance between two datasets.

        Parameters
        ----------
        d1, d2: list
            A list, containing a list describing the data set
        """
        if not isinstance(d1, list):
            raise TypeError('Data is not of allowed types')
        if not isinstance(d2, list):
            raise TypeError('Data is not of allowed types')
        
        # Distance main algorithm
        def single_distance(obsHeatmap,simHeatmap):
            resObsSorted = np.sort(np.reshape(obsHeatmap, (1, obsHeatmap.shape[0] * obsHeatmap.shape[1])))
            resSimSorted = np.sort(np.reshape(simHeatmap, (1, simHeatmap.shape[0] * simHeatmap.shape[1])))
            return sum(sum(pow(resObsSorted-resSimSorted, 2)))

        # Extract summary statistics from the dataset
        if (self.s1 is None or (np.asarray(self.data_set) == np.asarray(d1)).all() == False):
            self.s1 = self.statistics_calc.statistics(d1)
            self.data_set = d1
        s2 = self.statistics_calc.statistics(d2)

        ## EXTRACT MATRICES FROM LIST
        obsRep, simRep = len(self.s1), len(s2)

        # Observation
        resObservationHeatmap = []
        for repi in range(obsRep):
            # Extract fundamental information
            nCells = int(self.s1[repi,1])
            # Extract observation heat map
            repObservationHeatmap = np.zeros((nCells,nCells))
            for i in range(nCells):
                index = 2 + i*nCells
                repObservationHeatmap[i,:] = self.s1[repi,index:index+nCells]
            # Add to global quantities
            resObservationHeatmap.append(copy.deepcopy(repObservationHeatmap))
        
        # Simulation
        resSimulationHeatmap = []
        for repi in range(simRep):
            # Extract fundamental information
            nCells = int(s2[repi,1])
            # Extract simulation heat map
            repSimulationHeatmap = np.zeros((nCells,nCells))
            for i in range(nCells):
                index = 2 + i*nCells
                repSimulationHeatmap[i,:] = s2[repi,index:index+nCells]    
            # Add to global quantities
            resSimulationHeatmap.append(copy.deepcopy(repSimulationHeatmap))
        
        ## ACTUALLY COMPUTE DISTANCE
        # Compute distance using heatmap
        tempDistance = 0
        for simNumber in range(simRep):
            for obsNumber in range(obsRep):
                tempDistance += single_distance(resObservationHeatmap[obsNumber],resSimulationHeatmap[simNumber])
        return tempDistance/(nCells**2 * simRep * obsRep)

    def dist_max(self):
        return np.inf
        
class DistanceType3(Distance):
    """
    This is distance is more complex than the others, it is used to compare how different are the peaks both in size and location
    It has been self developed and aim at reproducing a peak similar to the one observed during the experiment

    The maximum value of the distance is np.inf.
    """

    def __init__(self, statistics):
        self.statistics_calc = statistics

        # Since the observations do always stay the same, we can save the
        #  summary statistics of them and not recalculate it each time
        self.s1 = None
        self.data_set = None
        
    def distance(self, d1, d2):
        """Calculates the distance between two datasets.

        Parameters
        ----------
        d1, d2: list
            A list, containing a list describing the data set
        """
        if not isinstance(d1, list):
            raise TypeError('Data is not of allowed types')
        if not isinstance(d2, list):
            raise TypeError('Data is not of allowed types')
        
        # Distance main algorithm
        def single_distance(obsHeatmap,simHeatmap):
            # Sort elements and get unique values
            resObsSorted = np.sort(np.reshape(obsHeatmap, (1, obsHeatmap.shape[0] * obsHeatmap.shape[1])))
            resSimSorted = np.sort(np.reshape(simHeatmap, (1, simHeatmap.shape[0] * simHeatmap.shape[1])))
            resObsElements = np.unique(resObsSorted)
            resSimElements = np.unique(resSimSorted)
            
            # Determine which matrix has most steps
            highIndex = np.array([resObsElements.shape[0],resSimElements.shape[0]]).argmax()
            if highIndex==0:    highArray, highElements, lowArray = obsHeatmap, resObsElements, simHeatmap
            else:               highArray, highElements, lowArray = simHeatmap, resSimElements, obsHeatmap
            
            distanceValue = 0
            for highValue in highElements:
                # Look for areas with same level in highArray and find closest height in lowArray
                highIndex, diff = [], 1E+10
                for i in range(highArray.shape[0]):
                    for j in range(highArray.shape[1]):
                        if highArray[i,j]==highValue:
                            highIndex.append(np.array([i,j]))
                        if abs(highValue-lowArray[i,j])<diff:
                            lowValue = lowArray[i,j]
                            diff = abs(highValue-lowArray[i,j])
                # Look for areas with similar level in lowArray
                lowIndex = []
                for i in range(lowArray.shape[0]):
                    for j in range(lowArray.shape[1]):
                        if lowArray[i,j]==lowValue:
                            lowIndex.append(np.array([i,j]))
                # Compute minimum distance between similar areas
                dist = 1E+10
                for i in range(len(lowIndex)):
                    for j in range(len(highIndex)):
                        diff = lowIndex[i]-highIndex[j]
                        if math.sqrt(diff[0]**2+diff[1]**2)<dist:
                            dist = math.sqrt(diff[0]**2+diff[1]**2)
                # Compute final distance
                distanceValue += abs(highValue-lowValue)*dist
            return distanceValue

        # Extract summary statistics from the dataset
        if (self.s1 is None or (np.asarray(self.data_set) == np.asarray(d1)).all() == False):
            self.s1 = self.statistics_calc.statistics(d1)
            self.data_set = d1
        s2 = self.statistics_calc.statistics(d2)

        ## EXTRACT MATRICES FROM LIST
        obsRep, simRep = len(self.s1), len(s2)

        # Observation
        resObservationHeatmap = []
        for repi in range(obsRep):
            # Extract fundamental information
            nCells = int(self.s1[repi,1])
            # Extract observation heat map
            repObservationHeatmap = np.zeros((nCells,nCells))
            for i in range(nCells):
                index = 2 + i*nCells
                repObservationHeatmap[i,:] = self.s1[repi,index:index+nCells]
            # Add to global quantities
            resObservationHeatmap.append(copy.deepcopy(repObservationHeatmap))
        
        # Simulation
        resSimulationHeatmap = []
        for repi in range(simRep):
            # Extract fundamental information
            nCells = int(s2[repi,1])
            # Extract simulation heat map
            repSimulationHeatmap = np.zeros((nCells,nCells))
            for i in range(nCells):
                index = 2 + i*nCells
                repSimulationHeatmap[i,:] = s2[repi,index:index+nCells]    
            # Add to global quantities
            resSimulationHeatmap.append(copy.deepcopy(repSimulationHeatmap))
        
        ## ACTUALLY COMPUTE DISTANCE
        # Compute distance using heatmap
        tempDistance = 0
        for simNumber in range(simRep):
            for obsNumber in range(obsRep):
                tempDistance += single_distance(resObservationHeatmap[obsNumber],resSimulationHeatmap[simNumber])
        return tempDistance/(nCells**2 * simRep * obsRep)

    def dist_max(self):
        return np.inf

class DistanceType4(Distance):
    """
    This distance compute the Euclidean distance between the positions of pedestrians.
    In this distance only 1 and 0 are considered.
    This distance was proposed by Antonietta Mira.

    The maximum value of the distance is np.inf.
    """

    def __init__(self, statistics):
        self.statistics_calc = statistics

        # Since the observations do always stay the same, we can save the
        #  summary statistics of them and not recalculate it each time
        self.s1 = None
        self.data_set = None

    def distance(self, d1, d2):
        """Calculates the distance between two datasets.

        Parameters
        ----------
        d1, d2: list
            A list, containing a list describing the data set
        """
        if not isinstance(d1, list):
            raise TypeError('Data is not of allowed types')
        if not isinstance(d2, list):
            raise TypeError('Data is not of allowed types')

        # Extract summary statistics from the dataset
        if (self.s1 is None or (np.asarray(self.data_set) == np.asarray(d1)).all() == False):
            self.s1 = self.statistics_calc.statistics(d1)
            self.data_set = d1
        s2 = self.statistics_calc.statistics(d2)

        ## EXTRACT MATRICES FROM LIST
        obsRep, simRep = len(self.s1), len(s2)
        
        # Observation
        resObservationTime, resObservationPos = [], []
        for repi in range(obsRep):
            # Extract fundamental information
            nTime = int(self.s1[repi,0])
            nCells = int(self.s1[repi,1])
            # Extract observation time
            repObservationTime = []
            for n in range(nTime):
                index = 2 + nCells**2 + n
                if n==0 or self.s1[repi,index]>0:
                    repObservationTime.append(copy.deepcopy(self.s1[repi,index]))
            # Extract observation positions
            repObservationPos = []
            for n in range(len(repObservationTime)):
                pedPos = np.zeros((nCells,nCells))
                for i in range(nCells):
                    index = 2 + nCells**2 + nTime + nCells**2 * n + i*nCells
                    pedPos[i,:] = self.s1[repi,index:index+nCells]
                repObservationPos.append(copy.deepcopy(pedPos))
            # Add to global quantities
            resObservationTime.append(copy.deepcopy(repObservationTime))
            resObservationPos.append(copy.deepcopy(repObservationPos))
        
        # Simulation
        resSimulationTime, resSimulationPos = [], []
        for repi in range(simRep):
            # Extract fundamental information
            nTime = int(s2[repi,0])
            nCells = int(s2[repi,1])
            # Extract simulation time
            repSimulationTime = []
            for n in range(nTime):
                index = 2 + nCells**2 + n
                if n==0 or s2[repi,index]>0:
                    repSimulationTime.append(copy.deepcopy(s2[repi,index]))
            # Extract simulation positions
            repSimulationPos = []
            for n in range(len(repSimulationTime)):
                pedPos = np.zeros((nCells,nCells))
                for i in range(nCells):
                    index = 2 + nCells**2 + nTime + nCells**2 * n + i*nCells
                    pedPos[i,:] = s2[repi,index:index+nCells]
                repSimulationPos.append(copy.deepcopy(pedPos))
            # Add to global quantities
            resSimulationTime.append(copy.deepcopy(repSimulationTime))
            resSimulationPos.append(copy.deepcopy(repSimulationPos))

        ## ACTUALLY COMPUTE DISTANCE
        # Compare times and get correspondings indexes
        indexLink = []
        for simNumber in range(len(resSimulationTime)):
            indexLinkLocal = []
            for obsNumber in range(len(resObservationTime)):
                tempLink, index = [], 0
                for i in range(len(resSimulationTime[simNumber])):
                    keepSearching = True
                    if index < len(resObservationTime[obsNumber]):
                        while keepSearching == True and resObservationTime[obsNumber][index] - \
                                resSimulationTime[simNumber][i] <= 1E-5:
                            index = index + 1
                            if index == len(resObservationTime[obsNumber]):
                                keepSearching = False
                        tempLink.append(copy.deepcopy(index - 1))
                indexLinkLocal.append(copy.deepcopy(tempLink))
            indexLink.append(copy.deepcopy(indexLinkLocal))

        # Compute distance using single positions
        totalSamples, diffMatrix = 0, np.zeros((nCells, nCells))
        for simNumber in range(simRep):
            for expNumber in range(obsRep):
                for n in range(len(indexLink[simNumber][expNumber])):
                    diffMatrix = diffMatrix + abs(resSimulationPos[simNumber][n] - resObservationPos[expNumber][
                        indexLink[simNumber][expNumber][n]])
                    totalSamples = totalSamples + 1
        return diffMatrix.sum() / (nCells ** 2 * totalSamples)

    def dist_max(self):
        return np.inf