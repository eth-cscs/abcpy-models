from abcpy.distances import Distance
import copy
import math
import numpy as np
from itertools import combinations
import statsmodels.formula.api as smf



class DistanceWeighted(Distance):
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
        self.D1 = DistanceType1()
        self.D2 = DistanceType2()
        self.D3 = DistanceType3()
        self.D4 = DistanceType4()
        self.weight = np.array([1, 1, 1, 1])

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

        Distances = np.array([self.D1.distance(d1, d2), self.D2.distance(d1, d2), self.D3.distance(d1, d2),
                     self.D4.distance(d1, d2)])
        Result = sum(self.weight * Distances)
        return Result

    def dist_max(self):
        Max = np.array([self.D1.dist_max(), self.D2.dist_max(), self.D3.dist_max(), self.D4.dist_max()])
        return sum(self.weight * Max)

    def update(self, new_all_parameters, new_all_data):
        param_distance, data_distance = [], []
        for element in list(combinations(range(len(new_all_parameters)), 2)):
            param_distance.append(sum(abs(np.concatenate(new_all_parameters[element[0]]) - np.concatenate(new_all_parameters[element[1]]))))
            data1, data2 = new_all_data[element[0]][0], new_all_data[element[1]][0]
            data_distance.append([self.D1.distance(data1, data2), self.D2.distance(data1, data2),
                                  self.D3.distance(data1, data2), self.D4.distance(data1, data2)])
        results = smf.OLS(np.array(param_distance), np.array(data_distance)).fit()
        self.weight = np.array(results.params)
        return self.weight


class DistanceType1(Distance):
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
                # Absolute distance between Heatmaps
                tempDistance += sum(sum(abs(resObservationHeatmap[obsNumber] - resSimulationHeatmap[simNumber]))) / pow(
                    nCell, 2)

        return tempDistance / (len(d1) * len(d2))

    def dist_max(self):
        return np.inf

class DistanceType2(Distance):
    """
    This distance is Euclidean distance between the sorted heatmaps
    The idea is that is not so important where the maximum value is but how big it is and how are the number distributed
    This is simply a guess and this distance was defined to check if such idea work or not

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
                # Sorted distance between heatmaps
                resObsSorted = np.sort(np.reshape(resObservationHeatmap[obsNumber], (1, resObservationHeatmap[obsNumber].shape[0] * resObservationHeatmap[obsNumber].shape[1])))
                simObsSorted = np.sort(np.reshape(resSimulationHeatmap[simNumber], (1, resSimulationHeatmap[simNumber].shape[0] * resSimulationHeatmap[simNumber].shape[1])))
                tempDistance += sum(sum(abs(resObsSorted - simObsSorted))) / pow(nCell, 2)
                # Sorted distance each time
                tempDistance += self.all_dist(resObsPos[obsNumber],resSimPos[simNumber])
                
        return tempDistance / (len(d1) * len(d2))

    def dist_max(self):
        return np.inf
    
    def all_dist(self, resObsPos, resSimPos):
        if len(resSimPos) > len(resObsPos):
            resObsPos = resObsPos.tolist() + [0 for i in range(len(resSimPos) - len(resObsPos))]
        else:
            resSimPos = resSimPos.tolist() + [0 for i in range(len(resObsPos) - len(resSimPos))]
        return sum(abs(resObsPos - resSimPos)) / len(resSimPos)
        
class DistanceType3(Distance):
    """
    This is distance is more complex than the others, it is used to compare how different are the peaks both in size and location
    It has been self developed and aim at reproducing a peak similar to the one observed during the experiment

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
            for obsNumber in range(len(d1)):
                resObsSorted = np.sort(np.reshape(resObservationHeatmap[obsNumber], (1, resObservationHeatmap[obsNumber].shape[0] * resObservationHeatmap[obsNumber].shape[1])))
                resSimSorted = np.sort(np.reshape(resSimulationHeatmap[simNumber], (1, resSimulationHeatmap[simNumber].shape[0] * resSimulationHeatmap[simNumber].shape[1])))
                resObsElements = np.unique(resObsSorted)
                resSimElements = np.unique(resSimSorted)
            
                # Determine which matrix has most steps
                highIndex = np.array([resObsElements.shape[0],resSimElements.shape[0]]).argmax()
                if highIndex==0:    highArray, highElements, lowArray = resObservationHeatmap[obsNumber], resObsElements, resSimulationHeatmap[simNumber]
                else:               highArray, highElements, lowArray = resSimulationHeatmap[simNumber], resSimElements, resObservationHeatmap[obsNumber]
                
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
                
                tempDistance += distanceValue / pow(nCell, 2)
                tempDistance += self.all_dist(resObsPos[obsNumber],resSimPos[simNumber])
                
        return tempDistance / (len(d1) * len(d2))

    def dist_max(self):
        return np.inf
    
    def all_dist(self, resObsPos, resSimPos):
        if len(resSimPos) > len(resObsPos):
            resObsPos = resObsPos.tolist() + [0 for i in range(len(resSimPos) - len(resObsPos))]
        else:
            resSimPos = resSimPos.tolist() + [0 for i in range(len(resObsPos) - len(resSimPos))]
        return sum(abs(resObsPos - resSimPos)) / len(resSimPos)

class DistanceType4(Distance):
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
                # Absolute Distance between pedpositions at each time
                tempDistance += self.all_dist(resObsPos[obsNumber], resSimPos[simNumber])

        return tempDistance / (len(d1) * len(d2))

    def dist_max(self):
        return np.inf

    def all_dist(self, resObsPos, resSimPos):
        if len(resSimPos) > len(resObsPos):
            resObsPos = resObsPos.tolist() + [0 for i in range(len(resSimPos) - len(resObsPos))]
        else:
            resSimPos = resSimPos.tolist() + [0 for i in range(len(resObsPos) - len(resSimPos))]
        return sum(abs(resObsPos - resSimPos)) / len(resSimPos)
