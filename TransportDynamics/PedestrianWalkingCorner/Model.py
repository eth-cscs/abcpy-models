import numpy as np
import math, csv, copy, sys
import random as rand

from abcpy.probabilisticmodels import ProbabilisticModel, Continuous, InputConnector, Discrete

class FFSimulator(ProbabilisticModel, Continuous):
    """
        Parameters
        ----------
        parameters: list
            Contains the probabilistic models and hyperparameters from which the model derives.

        """
    def __init__(self, parameters, name='FloorField'):

        # Parameter specifying the dimension of the return values of the distribution.
        #self.dimension = len(self.initial_state)
        # Other parameters kept fixed
        # Constants (related to experimental conditions, known a priori)
        self.maxPed = 16                    # number of pedestrians (whole group)
        self.vFree = 1.4                    # free walking speed (m/s)
        self.startArea = 16                 # starting area (m2)

        # Internal settings and options (determine details of the simulation)
        self.meshSize = 0.4             # mesh size (m)
        self.allowStop = 0              # allow stopping (if 0 stopping not allowed)
        self.motionLogic = 0            # consider priorities in deciding movements (right vs. left,...)
        self.maxInflowIter = 100        # iteration limits for inflow (avoid being stuck in a while loop)
        self.maxSimulationTime = 30     # time limit for simulation (s) (in case simulation get stuck due to odd parameter's choice)

        # We expect input of type parameters = [theta1, theta2, n_timestep]
        if not isinstance(parameters, list):
            raise TypeError('Input of FloorField model is of type list')

        if len(parameters) != 7:
            raise RuntimeError('Input list must be of length 7, containing [bufferRatio, bufferAngle, kW, kS, kD, decay, diffusion].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)


    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 7:
            raise ValueError('Number of parameters of FloorField model must be 7.')

        # Check whether input is from correct domain
        bufferRatio = input_values[0]       # length ratio of buffer region                         # PARAMETER
        bufferAngle = input_values[1]       # deflection angle in the buffer region (deg) [45,90)   # PARAMETER
        kW = input_values[2]                # wall floor field weighting factor                     # PARAMETER
        kS = input_values[3]                # static floor field weighting factor                   # parameter
        kD = input_values[4]                # dynamic floor field weighting factor                  # parameter
        decay = input_values[5]             # decay of the dynamic FF                               # parameter
        diffusion = input_values[6]         # diffusion of the dynamic FF                           # parameter

        if not isinstance(bufferRatio, (float, np.float64, np.float32, np.float16)) or bufferRatio < 0:
            print('bufferRatio is not of correct type or out of range')
            return False

        if not isinstance(bufferAngle, (int, np.int64, np.int32, np.int16)) or bufferAngle < 45 or bufferAngle > 90:
            print('bufferAngle is not of correct type or out of range')
            return False

        if not isinstance(kW, (float, np.float64, np.float32, np.float16)) or kW < 0:
            print('kW is not of correct type or out of range')
            return False

        if not isinstance(kS, (float, np.float64, np.float32, np.float16)) or kS < 0:
            print('kS is not of correct type or out of range')
            return False

        if not isinstance(kD, (float, np.float64, np.float32, np.float16)) or kD < 0:
            print('kD is not of correct type or out of range')
            return False

        if not isinstance(decay, (float, np.float64, np.float32, np.float16)) or decay < 0 or diffusion > 1:
            print('decay is not of correct type or out of range')
            return False

        if not isinstance(diffusion, (float, np.float64, np.float32, np.float16)) or diffusion < 0 or diffusion > 1:
            print('diffusion is not of correct type or out of range')
            return False

        return True


    def _check_output(self, values):
        if not isinstance(values[0], np.ndarray):
            raise ValueError('Output of the normal distribution is always a number.')
        return True

    def get_output_dimension(self):
        return 1

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        # Extract the input parameters
        bufferRatio = input_values[0]       # length ratio of buffer region                         # PARAMETER
        bufferAngle = input_values[1]       # deflection angle in the buffer region (deg) [45,90)   # PARAMETER
        kW = input_values[2]                # wall floor field weighting factor                     # PARAMETER
        kS = input_values[3]                # static floor field weighting factor                   # parameter
        kD = input_values[4]                # dynamic floor field weighting factor                  # parameter
        decay = input_values[5]             # decay of the dynamic FF                               # parameter
        diffusion = input_values[6]         # diffusion of the dynamic FF                           # parameter

        # Do the actual forward simulation
        vector_of_k_samples = self.floorfield(bufferRatio, bufferAngle, kW, kS, kD, decay, diffusion, k)
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result


    def floorfield(self, bufferRatio, bufferAngle, kW, kS, kD, decay, diffusion, simRep):
        # Generate n_simulated data
        result = []

        ## CREATE STATIC AND WALL FLOOR FIELDS (both depends on parameters chosen)
        # In floorMap: 0 is wall or obstacles       -1 is entrance       2 is a special region not considered in results
        #              1 is walkable area           -2 is exit             (it is used to directly push people to the exit)

        # Read floor map and determine section dimension and corridor width
        floorMap = np.array(list(csv.reader(open("floor_map/Lcorner.csv"), delimiter=","))).astype("float")
        (iSize, jSize) = floorMap.shape
        entranceLocation = [math.inf, 0]
        for i in range(iSize):
            if floorMap[i][jSize - 1] == -1:
                entranceLocation[0] = min([i, entranceLocation[0]]);
                entranceLocation[1] = max([i, entranceLocation[1]])
        corrWidth = entranceLocation[1] - entranceLocation[0] + 1

        # Determine location of turning point
        maxDistCenter, maxDistCorner = 0, 0
        iCenter, jCenter, iCorner, jCorner = 0, 0, 0, 0
        for i in range(iSize - 1, 0, -1):
            for j in range(jSize - 1, 0, -1):
                dist = (iSize - i) ** 2 + (jSize - j) ** 2
                if dist > maxDistCenter and floorMap[i][j] == 0:        maxDistCenter, iCenter, jCenter = dist, i, j
                if dist > maxDistCorner and floorMap[i][j] == 1:        maxDistCorner, iCorner, jCorner = dist, i, j

        # Create floor map for final approach (last straight part)
        staticFF = np.zeros((iSize, jSize))
        for i in range(iSize - 1, 0, -1):
            if i >= iCenter:      staticFF[i, 1:corrWidth + 1] = iSize - i - 1

        # Create radial map for turning point (central part with curve)
        startValue = np.amax(staticFF)
        for i in range(iSize):
            for j in range(jSize):
                if i < iCenter and j < jCenter and i > 0 and j > 0:
                    angle = math.atan2((i - iCenter) * self.meshSize, (jCenter - j) * self.meshSize)
                    staticFF[i][j] = abs(angle * (corrWidth / (math.pi / 2))) + startValue

        # Create deflection before turning (according to Li et al., Physica A 432 (2015) 337–353)
        bufferCells, slope = round(bufferRatio * corrWidth), math.tan(abs(bufferAngle - 90) * (math.pi / 180))

        def getDistance(slope, i, j):
            return abs(iCenter - slope * jCenter + slope * j - i) / math.sqrt(1 + slope ** 2)

        minDist, maxDist = getDistance(slope, iCenter - 1, jCenter), getDistance(slope, iCenter - corrWidth,
                                                                                 jCenter + bufferCells - 1)
        for i in range(iSize):
            for j in range(jSize):
                if i < iCenter and i > 0 and j >= jCenter and j < jCenter + bufferCells:
                    bufferValue = getDistance(slope, i, j)
                    staticFF[i, j] = (bufferValue - minDist) * (
                    bufferCells / (maxDist - minDist)) + jSize - jCenter + corrWidth - 2

        # Create floor map for initial approach (initial straight section)
        startValue = np.amax(staticFF) - jCenter - bufferCells + 1
        for j in range(jCenter, jSize):
            if j >= jCenter + bufferCells:      staticFF[1:corrWidth + 1, j] = startValue + j

        # Normalize to have only integer values
        allData = staticFF.flatten()
        setValues = sorted(set(allData))
        for i in range(iSize):
            for j in range(jSize):
                staticFF[i, j] = setValues.index(staticFF[i, j])

        # Compute wall floor field (Nishinari et al., IEICE Transactions on information and systems 87.3, 726-732)
        wallFF = np.zeros((iSize, jSize))
        for i in range(iSize):
            for j in range(jSize):
                if floorMap[i][j] != 0:
                    minDist = math.inf
                    for m in range(iSize):
                        for n in range(jSize):
                            if floorMap[m][n] == 0:
                                dist = math.sqrt((m - i) ** 2 + (n - j) ** 2)
                                if dist < minDist:    minDist = dist
                    wallFF[i][j] = minDist

        minExp, maxExp = math.log(sys.float_info.min), math.log(sys.float_info.max)

        def compProb(sFF, dFF, wFF):
            exponent = -kS * sFF + kD * dFF + kW * wFF
            if exponent > minExp and exponent < maxExp:
                return math.exp(-kS * sFF + kD * dFF + kW * wFF)
            elif exponent < minExp:
                return math.exp(minExp)
            elif exponent > maxExp:
                return math.exp(maxExp)

        compProb = np.vectorize(compProb)

        # Main simulation loop  (including repetitions)
        for repi in range(simRep):
            
            # Initialize variables (relative to single execution)
            repSimulationPos, repSimulationTime = [], []
            pedPos, repSimulationHeatmap = np.zeros((iSize, jSize)), np.zeros((iSize-3,jSize-3))
            dynamicFF = np.zeros((iSize, jSize))
            pedList, pedDir = [], []
            numPed, inTime, simTime, tZero = 0, 0, 0, 0

            # Internal simulation loop (single run)
            while ((numPed < self.maxPed) or (len(pedList) > 0)) and (simTime < self.maxSimulationTime):
                # Inject pedestrians if required
                if numPed < self.maxPed:
                    warningCounter = 0
                    while (numPed < self.vFree * corrWidth * self.meshSize * (self.maxPed / self.startArea) * (simTime - inTime)):
                        m, n = rand.randint(entranceLocation[0], entranceLocation[1]), jSize - 1
                        if pedPos[m][n] == 0:
                            pedPos[m][n] = 1
                            numPed = numPed + 1
                            pedList.append([m, n])
                            pedDir.append([0, -1])
                        if numPed == 1:
                            inTime = simTime
                        warningCounter = warningCounter + 1
                        if warningCounter > self.maxInflowIter:      numPed = numPed + 1

                # Remove pedestrians if required
                remPed = []
                for pedID in range(len(pedList)):
                    if abs(floorMap[pedList[pedID][0]][pedList[pedID][1]]) == 2:      remPed.append(pedID)
                remPed.sort(reverse=True)
                for pedID in range(len(remPed)):
                    del pedList[remPed[pedID]]
                    del pedDir[remPed[pedID]]

                # Compute transition probabilities and reserve potential positions
                pedRes = [[[[] for j in range(jSize)] for i in range(iSize)] for n in range(2)]
                for pedID in range(len(pedList)):

                    # Extract local floor fields and prepare data for transition probability calculation
                    I, J = pedList[pedID][0], pedList[pedID][1]
                    localMap = [None] * 3
                    localMap[0] = compProb(staticFF[I - 1:I + 2, J - 1:J + 2], dynamicFF[I - 1:I + 2, J - 1:J + 2],
                                           wallFF[I - 1:I + 2, J - 1:J + 2])
                    localMap[1], localMap[2] = pedPos[I - 1:I + 2, J - 1:J + 2], floorMap[I - 1:I + 2, J - 1:J + 2]
                    for n in range(len(localMap)):
                        if J == jSize - 1:
                            localMap[n] = np.c_[localMap[n], np.zeros(3)]
                        elif J == 0:
                            localMap[n] = np.c_[np.zeros(3), localMap[n]]
                        elif I == iSize - 1:
                            localMap[n] = np.r_[localMap[n], [np.zeros(3)]]
                        elif I == 0:
                            localMap[n] = np.r_[[np.zeros(3)], localMap[n]]

                    # Exclude corners (Von Neumann neighborhood), back step, stop (if not allowed) and inacessible locations
                    localMap[0][0][0], localMap[0][2][0], localMap[0][0][2], localMap[0][2][2] = 0, 0, 0, 0
                    dI, dJ = pedDir[pedID][0], pedDir[pedID][1]
                    localMap[0][1 - dI][1 - dJ] = 0  # remove back stepping
                    if self.allowStop == 0:        localMap[0][1][1] = 0  # no stopping allows (depending on options)
                    for i in range(localMap[0].shape[0]):
                        for j in range(localMap[0].shape[1]):
                            if not (i == 1 and j == 1) and localMap[1][i][j] == 1:     localMap[0][i][
                                j] = 0  # exclude occupied positions
                            if localMap[2][i][j] == 0:                            localMap[0][i][
                                j] = 0  # exclude inacessible locations
                    totalValue = sum(sum(localMap[0]))
                    if totalValue > 0:    transProb = localMap[0] / totalValue  # compute transition probability

                    # Determine maximum values and relative index (also for tied cases)
                    maxValue, maxIndex = transProb.max(), []
                    if totalValue > 0:
                        for i in range(transProb.shape[0]):
                            for j in range(transProb.shape[1]):
                                if transProb[i][j] == maxValue:
                                    maxIndex.append([i, j])
                        if len(maxIndex) == 1:
                            m, n = maxIndex[0]  # only one maximum
                        else:
                            if self.motionLogic == 0:
                                m, n = maxIndex[rand.randint(0, len(maxIndex) - 1)]
                            else:
                                if len(maxIndex) == 4:
                                    m, n = 1 + dI, 1 + dJ  # all tied, keep the same walking direction
                                elif len(maxIndex) == 2:
                                    if [1 + dI, 1 + dJ] in maxIndex:
                                        m, n = 1 + dI, 1 + dJ  # keep walking direction
                                    elif maxIndex[0] == [1, 1]:
                                        m, n = maxIndex[1]  # do not stop (choose other direction)
                                    elif maxIndex[1] == [1, 1]:
                                        m, n = maxIndex[0]  # do not stop (choose other direction)
                                    else:
                                        m, n = 1 + dJ, 1 - dI  # turn right
                                elif len(maxIndex) == 3:
                                    if [1 + dI, 1 + dJ] in maxIndex:
                                        m, n = 1 + dI, 1 + dJ  # keep walking direction
                                    else:
                                        m, n = 1 + dJ, 1 - dI  # turn right
                    else:
                        m, n = 1, 1  # special case where transition probability is 0 everywhere

                    # Reserve position
                    pedRes[0][I + (m - 1)][J + (n - 1)].append(pedID)
                    if floorMap[I + (m - 1)][J + (n - 1)] == 2:
                        pedRes[1][I + (m - 1)][J + (n - 1)].append([dI, dJ])
                    else:
                        pedRes[1][I + (m - 1)][J + (n - 1)].append([m - 1, n - 1])

                # Resolve conflicts (in areas with floorMap==2 conflics are possible to avoid congestion near exit)
                for i in range(iSize):
                    for j in range(jSize - 2):
                        nReserved = len(pedRes[0][i][j])
                        if nReserved > 1:
                            chosenIndex, remPed = rand.randint(0, nReserved - 1), []
                            for n in range(nReserved):
                                if n != chosenIndex:      remPed.append(pedRes[0][i][j][n])
                            pedRes[0][i][j] = [pedRes[0][i][j][chosenIndex]]
                            pedRes[1][i][j] = [pedRes[1][i][j][chosenIndex]]
                            for n in range(len(remPed)):
                                pedID = remPed[n]
                                I, J = pedList[pedID][0], pedList[pedID][1]
                                pedRes[0][I][J] = [pedID]
                                pedRes[1][I][J] = [[0, 0]]

                # Move pedestrian and update position, direction and dynamic floor field
                pedPos = np.zeros((iSize, jSize))
                pedMoved = np.zeros((iSize,jSize))
                for i in range(iSize):
                    for j in range(jSize):
                        for n in range(len(pedRes[0][i][j])):
                            pedID = pedRes[0][i][j][n]
                            I, J = pedList[pedID][0], pedList[pedID][1]
                            dynamicFF[I][J] = dynamicFF[I][J] + 1
                            pedList[pedID] = [i, j]
                            pedDir[pedID] = pedRes[1][i][j][n]
                        if len(pedRes[0][i][j]) > 0:
                            pedPos[i][j] = 1
                            if abs(pedDir[pedID][0])+abs(pedDir[pedID][1])>0:
                                pedMoved[i][j] = 1

                # Diffuse and decay dynamic floor field
                for i in range(iSize):
                    for j in range(jSize):
                        if floorMap[i][j] != 0:
                            if dynamicFF[i][j] > 0 and rand.uniform(0, 1) < decay:
                                dynamicFF[i][j] = dynamicFF[i][j] - 1
                            if dynamicFF[i][j] > 0 and rand.uniform(0, 1) < diffusion:
                                dynamicFF[i][j] = dynamicFF[i][j] - 1
                                dI = int(np.sign(rand.uniform(-1, 1)))
                                dJ = int(np.sign(rand.uniform(-1, 1)))
                                if i + dI < iSize and j + dJ < jSize:
                                    dynamicFF[i + dI][j + dJ] = dynamicFF[i + dI][j + dJ] + 1

                # Store data for current iteration
                simTime = simTime + (self.meshSize / self.vFree)
                if sum(sum(pedPos[0:-3, 0:-1])) > 0:
                    if tZero == 0:        tZero = simTime
                    repSimulationPos.append(copy.deepcopy(pedPos[0:-3, 0:-3]))
                    repSimulationTime.append(copy.deepcopy(simTime - tZero))
                    repSimulationHeatmap = repSimulationHeatmap+pedMoved[0:-3,0:-3] * (self.meshSize / self.vFree)
            
            #We put repSimulationPos, repSimulationTime and repSimulationHeatmap in a list. It's of size (1, 2+p^2+nt+(p^2)*nt)
            #Notice we would expect all datasets simulated or observed should be stored like this and while computing distance
            #we will use this convention. - This is done as ABCpy expects a dataset as a nparray.
            nTime = math.ceil(self.maxSimulationTime / ((self.meshSize / self.vFree)))
            nCells = np.array(repSimulationHeatmap).shape[0]
            totElements = 2 + nCells**2 + nTime + (nCells**2) * nTime
            # Addiing primary information
            repResult = list([nTime, nCells])
            # Adding the SimulationHeatmap
            repResult += list(np.array(repSimulationHeatmap).reshape(-1,))
            # Adding the SimulationTime
            repResult += list(repSimulationTime)
            for inadd in range(len(repResult),2 + nCells**2 + nTime):
                repResult += list(np.array([0]))
            # Adding the SimulationPos
            for indadd in range(len(repSimulationPos)):
                repResult += list(np.array(repSimulationPos[indadd]).reshape(-1,))
            # Fill remaining elements
            for indadd in range(len(repResult),totElements):
                repResult += list(np.array([0]))

            # Append the result one simulation to the final simulation result
            result.append(copy.deepcopy(repResult))
      
        return result


class DiscreteUniform(Discrete, ProbabilisticModel):
    def __init__(self, parameters, name='DiscreteUniform'):
        """This class implements a probabilistic model following a Discrete Uniform distribution.

        Parameters
        ----------
        parameters: list
             A list containing two entries, the upper and lower bound of the range.

        name: string
            The name that should be given to the probabilistic model in the journal file.
        """

        if not isinstance(parameters, list):
            raise TypeError('Input for Discrete Uniform has to be of type list.')
        if len(parameters) != 2:
            raise ValueError('Input for Discrete Uniform has to be of length 2.')

        self._dimension = 1
        input_parameters = InputConnector.from_list(parameters)
        super(DiscreteUniform, self).__init__(input_parameters, name)
        self.visited = False

    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 2:
            raise ValueError('Number of parameters of FloorField model must be 2.')

        # Check whether input is from correct domain
        lowerbound = input_values[0]  # Lower bound
        upperbound = input_values[1]  # Upper bound

        if not isinstance(lowerbound, (int, np.int64, np.int32, np.int16)) or not isinstance(upperbound, (int, np.int64, np.int32, np.int16)) or lowerbound > upperbound:
            print('Parameters are not of correct type or out of range')
            return False

        return True

    def _check_output(self, parameters):
        """
        Checks parameter values given as fixed values. Returns False iff it is not an integer
        """
        if not isinstance(parameters[0], (int, np.int32, np.int64)):
            return False
        return True

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        """
        Samples from the Discrete Uniform distribution associated with the probabilistic model.

        Parameters
        ----------
        input_values: list
            List of input parameters, in the same order as specified in the InputConnector passed to the init function
        k: integer
            The number of samples to be drawn.
        rng: random number generator
            The random number generator to be used.

        Returns
        -------
        list: [np.ndarray]
            A list containing the sampled values as np-array.
        """

        result = np.array(rng.randint(input_values[0], input_values[1], size=k, dtype=np.int64))
        return [np.array([x]) for x in result]

    def get_output_dimension(self):
        return self._dimension

    def pmf(self, input_values, x):
        """Evaluates the probability mass function at point x.

        Parameters
        ----------
        input_values: list
            List of input parameters, in the same order as specified in the InputConnector passed to the init function
        x: float
            The point at which the pmf should be evaluated.

        Returns
        -------
        float:
            The pmf evaluated at point x.
        """
        lowerbound, upperbound = input_values[0], input_values[1]
        pmf = 1. / (upperbound - lowerbound + 1)
        self.calculated_pmf = pmf
        return pmf

