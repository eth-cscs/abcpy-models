from abcpy.distances import Distance
import copy
import numpy as np

class DistanceType1(Distance):
    """
    This class implements the Euclidean distance between two vectors.

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
        if (self.s1 is None or self.data_set != d1):
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
        tempDistance = np.zeros((1,simRep))
        for simNumber in range(simRep):
            diffMatrix = np.zeros((nCells,nCells))
            for obsNumber in range(obsRep):
                diffMatrix = diffMatrix+abs(sum(sum(pow(resSimulationHeatmap[simNumber]-resObservationHeatmap[obsNumber],2))))
            tempDistance[0,simNumber] = diffMatrix.sum()/(nCells**2 * obsRep)
        return np.mean(tempDistance)

    def dist_max(self):
        return np.inf


class DistanceType2(Distance):
    """
    This class implements the Euclidean distance between two vectors.

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
        if (self.s1 is None or self.data_set != d1):
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