import numpy as np
import math, csv, copy, sys
import random as rand
from scipy.stats import beta

from abcpy.probabilisticmodels import ProbabilisticModel, Continuous, InputConnector, Discrete

class FFSimulator1(ProbabilisticModel, Continuous):
    """
        Simulation using a simple approach proposed by Katsuhiro Nishinari
        
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
        self.totalPed = 42                  # number of pedestrians (whole group)
        self.vFree = 1.4                    # free walking speed (m/s)
        self.startSurface = 16              # starting area (m2)

        # Internal settings and options (determine details of the simulation)
        self.motionType = 2             # rule used to determine motion (1 = deterministic approach, 2 = probabilistic approach)
        self.entranceType = 1           # method used to introduce pedestrians (1 = set flow, 2 = use waiting area)
        self.meshSize = 0.4             # mesh size (m)
        self.allowStop = 0              # allow stopping (if 0 stopping not allowed)
        self.maxInflowIter = 100        # iteration limits for inflow (avoid being stuck in a while loop)
        self.maxSimulationTime = 60     # time limit for simulation (s) (in case simulation get stuck due to odd parameter's choice)

        # We expect input of type parameters = [theta1, theta2, n_timestep]
        if not isinstance(parameters, list):
            raise TypeError('Input of FloorField model is of type list')

        if len(parameters) != 5:
            raise RuntimeError('Input list must be of length 5, containing [kS, kD, kW, decay, diffusion].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)


    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 5:
            raise ValueError('Number of parameters of FloorField model must be 5.')

        # Check whether input is from correct domain
        kS = input_values[0]                # static floor field weighting factor
        kD = input_values[1]                # dynamic floor field weighting factor
        kW = input_values[2]                # wall floor field weighting factor
        decay = input_values[3]             # decay of the dynamic FF
        diffusion = input_values[4]         # diffusion of the dynamic FF
        
        if not isinstance(kS, (float, np.float64, np.float32, np.float16)) or kS < 0:
            print('kS is not of correct type or out of range')
            return False
        
        if not isinstance(kD, (float, np.float64, np.float32, np.float16)) or kD < 0:
            print('kD is not of correct type or out of range')
            return False
        
        if not isinstance(kW, (float, np.float64, np.float32, np.float16)) or kW < 0:
            print('kW is not of correct type or out of range')
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
        kS = input_values[0]                # static floor field weighting factor
        kD = input_values[1]                # dynamic floor field weighting factor
        kW = input_values[2]                # wall floor field weighting factor
        decay = input_values[3]             # decay of the dynamic FF
        diffusion = input_values[4]         # diffusion of the dynamic FF

        # Do the actual forward simulation
        vector_of_k_samples = self.floorfield(kS, kD, kW, decay, diffusion, k)
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result


    def floorfield(self, kS, kD, kW, decay, diffusion, simRep):
        # Generate n_simulated data
        result = []

        ## CREATE STATIC AND WALL FLOOR FIELDS (both depends on parameters chosen)
        # In floorMap: 0 is wall or obstacles       -1 is entrance       2 is a special region not considered in results
        #              1 is walkable area           -2 is exit             (it is used to directly push people to the exit)

        # Read floor map and determine section dimension and corridor width
        if self.entranceType==1:
            floorMap = np.array(list(csv.reader(open("floor_map/Lcorner.csv"), delimiter=","))).astype("float")
            (iSize, jSize) = floorMap.shape
            entranceLocation = [math.inf, 0]
            for i in range(iSize):
                if floorMap[i][jSize - 1] == -1:
                    entranceLocation[0] = min([i, entranceLocation[0]]);
                    entranceLocation[1] = max([i, entranceLocation[1]])
            corrWidth = entranceLocation[1] - entranceLocation[0] + 1
            
        if self.entranceType==2:
            floorMap = np.array(list(csv.reader(open("floor_map/Lcorner_long.csv"), delimiter=","))).astype("float")
            (iSize,jSize) = floorMap.shape
            corrWidth = int(abs(sum(floorMap[:,-1])))
            startArea = [[iSize,0],[jSize,0]]
            for i in range(iSize):
                for j in range(jSize):
                    if floorMap[i,j]==-1:
                        if startArea[0][0]>i:  startArea[0][0] = i
                        if startArea[0][1]<i:  startArea[0][1] = i
                        if startArea[1][0]>j:  startArea[1][0] = j
                        if startArea[1][1]<j:  startArea[1][1] = j

        # Determine location of turning point
        maxDistCenter, iCenter, jCenter = 0, 0, 0
        for i in range(iSize - 1, 0, -1):
            for j in range(jSize - 1, 0, -1):
                dist = (iSize - i) ** 2 + (jSize - j) ** 2
                if dist > maxDistCenter and floorMap[i][j] == 0:        maxDistCenter, iCenter, jCenter = dist, i, j
        
        # Static Floor Field using a simple approach by Nishinari         
        # Create floor map for last straight section
        staticFF = np.zeros((iSize,jSize))
        for i in range(iSize-1,0,-1):
            if i>=1:        staticFF[i,1:corrWidth+1] = iSize-i-1     
                
        # Create floor map for first straight section
        for j in range(1,jSize):
            if staticFF[1,j]==0:    staticFF[1:corrWidth+1,j] = staticFF[iCenter-1,jCenter-1]+(j-jCenter)
                
        # Create floor field for corner section
        for i in range(corrWidth):
            for j in range(corrWidth):
                if i>=j:
                    staticFF[iCenter-i-1,jCenter+j] += i-j

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

        minExp, maxExp = math.log(sys.float_info.min*5), math.log(sys.float_info.max/5)

        def compProb(sFFcenter, sFF, dFF, wFF):
            exponent = -kS * (sFF - sFFcenter) + kD * dFF + kW * wFF
            if exponent > minExp and exponent < maxExp:
                return math.exp(-kS * (sFF - sFFcenter) + kD * dFF + kW * wFF)
            elif exponent < minExp:
                return math.exp(minExp)
            elif exponent > maxExp:
                return math.exp(maxExp)

        compProb = np.vectorize(compProb)

        # Main simulation loop  (including repetitions)
        for repi in range(simRep):
            
            # Initialize variables (relative to single execution)
            repSimulationPos, repSimulationTime = [], []
            pedPos, repSimulationHeatmap = np.zeros((iSize, jSize)), np.zeros((iSize-3,iSize-3))
            dynamicFF = np.zeros((iSize, jSize))
            pedList, pedDir = [], []
            numPed, inTime, simTime, tZero = 0, 0, 0, 0
            
            # Place pedestrians in starting area
            if self.entranceType==2:
                while (numPed < self.totalPed):
                    m, n = rand.randint(startArea[0][0],startArea[0][1]), rand.randint(startArea[1][0],startArea[1][1])
                    if pedPos[m,n]==0:
                        pedPos[m][n] = 1
                        numPed = numPed + 1
                        pedList.append([m,n])
                        pedDir.append([0,-1])

            # Internal simulation loop (single run)
            maxTime = self.maxSimulationTime - (self.meshSize/self.vFree)
            while ((numPed < self.totalPed) or (len(pedList) > 0)) and (simTime-inTime < maxTime):
                
                # Inject pedestrians if required (entrance by inflow value)
                if self.entranceType==1:
                    if numPed < self.totalPed:
                        warningCounter = 0
                        while (numPed < self.vFree * corrWidth * self.meshSize * (self.totalPed / self.startSurface) * (simTime - inTime)):
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

                 # Define starting time for starting area setup
                if self.entranceType==2:
                    if sum(pedPos[:,iSize])>1 and inTime==0:
                        inTime = simTime                
                
                # Remove pedestrians if required
                remPed = []
                for pedID in range(len(pedList)):
                    if abs(floorMap[pedList[pedID][0]][pedList[pedID][1]]) == 2:      remPed.append(pedID)
                remPed.sort(reverse=True)
                for pedID in range(len(remPed)):
                    del pedList[remPed[pedID]]
                    del pedDir[remPed[pedID]]
                    
                # Diffuse and decay dynamic floor field
                for i in range(iSize):
                    for j in range(jSize):
                        if floorMap[i][j] != 0:
                            # Decay
                            localDynamicFF = int(dynamicFF[i][j])
                            if dynamicFF[i][j] > 0:
                                for n in range(0,localDynamicFF):
                                    if rand.uniform(0,1) < decay:
                                        dynamicFF[i][j] = dynamicFF[i][j] - 1
                            # Diffusion
                            if dynamicFF[i][j] > 0 and rand.uniform(0,1) < diffusion:
                                dynamicFF[i][j] = dynamicFF[i][j] - 1
                                index = rand.randint(0,1)
                                if index==0:
                                    dI, dJ = int(np.sign(rand.uniform(-1,1))), 0
                                elif index==1:
                                    dI, dJ = 0, int(np.sign(rand.uniform(-1,1)))
                                if i+dI < iSize and j+dJ < jSize and floorMap[i+dI][j+dJ]!=0:
                                    dynamicFF[i+dI][j+dJ] = dynamicFF[i+dI][j+dJ] + 1    

                # Compute transition probabilities and reserve potential positions
                pedRes = [[[[] for j in range(jSize)] for i in range(iSize)] for n in range(2)]
                for pedID in range(len(pedList)):
                    
                    # Extract local floor fields and prepare data for transition probability calculation
                    I, dI, J, dJ = pedList[pedID][0], pedDir[pedID][0], pedList[pedID][1], pedDir[pedID][1]
                    localMap = [None]*3
                    localMap[0] = compProb(staticFF[I,J], staticFF[I-1:I+2,J-1:J+2], dynamicFF[I-1:I+2,J-1:J+2], wallFF[I-1:I+2,J-1:J+2])
                    localMap[1], localMap[2] = pedPos[I-1:I+2,J-1:J+2], floorMap[I-1:I+2,J-1:J+2]
                    if (dI!=0 or dJ!=0) and I-dI>=0 and I-dI<iSize and J-dJ>=0 and J-dJ<jSize:
                        localValue = compProb(staticFF[I,J], staticFF[I-dI,J-dJ], dynamicFF[I-dI,J-dJ]-1, wallFF[I-dI,J-dJ])
                        localMap[0][1-dI][1-dJ] = localValue.item(0)
                    for n in range(len(localMap)):
                        if J==jSize-1:      localMap[n] = np.c_[localMap[n], np.zeros(3)]
                        elif J==0:          localMap[n] = np.c_[np.zeros(3), localMap[n]]
                        elif I==iSize-1:    localMap[n] = np.r_[localMap[n], [np.zeros(3)]]
                        elif I==0:          localMap[n] = np.r_[[np.zeros(3)], localMap[n]]
                        
                    # Exclude corners (Von Neumann neighborhood), back step, stop (if not allowed) and inacessible locations
                    localMap[0][0][0], localMap[0][2][0], localMap[0][0][2], localMap[0][2][2] = 0, 0, 0, 0
                    
                    localMap[0][1-dI][1-dJ] = 0                             # remove back stepping
                    if self.allowStop==0:        localMap[0][1][1] = 0      # no stopping allowed (depending on options)
                    for i in range(localMap[0].shape[0]):
                        for j in range(localMap[0].shape[1]):
                            if not(i==1 and j==1) and localMap[1][i][j]==1:     localMap[0][i][j] = 0       # exclude occupied positions
                            if localMap[2][i][j]==0:                            localMap[0][i][j] = 0       # exclude inacessible locations   
                    totalValue = sum(sum(localMap[0]))
                    
                    # Determine cell to reserve for next step
                    if self.motionType==1:
                        # Determine maximum values and relative index (also for tied cases)
                        if totalValue>0:    transProb = localMap[0]/totalValue          # compute transition probability
                        maxValue, maxIndex = transProb.max(), []
                        if totalValue>0:
                            for i in range(transProb.shape[0]):
                                for j in range(transProb.shape[1]):
                                    if transProb[i][j]==maxValue:
                                        maxIndex.append([i,j])
                            if len(maxIndex)==1:        m, n = maxIndex[0]      # only one maximum
                            else:                       m, n = maxIndex[rand.randint(0,len(maxIndex)-1)]
                        else:   m, n = 1, 1         # special case where transition probability is 0 everywhere
                    elif self.motionType==2:
                        # Determine next cell where do move stochastically based on previous results
                        p0, p1, p2, p3, p4 = localMap[0][1][1], localMap[0][1][0], localMap[0][1][2], localMap[0][0][1], localMap[0][2][1]
                        if totalValue>0:
                            randValue = rand.uniform(0,1)
                            if randValue >= (p0/totalValue) and randValue < ((p0+p1)/totalValue):
                                m, n = 1, 0
                            elif randValue >= ((p0+p1)/totalValue) and randValue < ((p0+p1+p2)/totalValue):
                                m, n = 1, 2
                            elif randValue >= ((p0+p1+p2)/totalValue) and randValue < ((p0+p1+p2+p3)/totalValue):
                                m, n = 0, 1
                            elif randValue >= ((p0+p1+p2+p3)/totalValue) and randValue < ((p0+p1+p2+p3+p4)/totalValue):
                                m, n = 2, 1
                            else:
                                m, n = 1, 1
                        else:
                            m, n = 1, 1
                        
                    # Reserve position
                    if I+(m-1)>=0 and I+(m-1)<iSize and J+(n-1)>=0 and J+(n-1)<jSize:
                        pedRes[0][I+(m-1)][J+(n-1)].append(pedID)
                        if floorMap[I+(m-1)][J+(n-1)]==2:   pedRes[1][I+(m-1)][J+(n-1)].append([dI,dJ])
                        else:                               pedRes[1][I+(m-1)][J+(n-1)].append([m-1,n-1])
                    
                # Resolve conflicts (in areas with floorMap==2, conflics are possible to avoid congestion near exit)
                for i in range(iSize):
                    for j in range(jSize-2):
                        nReserved = len(pedRes[0][i][j])
                        if nReserved>1:
                            # Determine winner candidate
                            chosenIndex, remPed = rand.randint(0,nReserved-1), []
                            # Losers are listed to move them back later
                            for n in range(nReserved):
                                if n!=chosenIndex:      remPed.append(pedRes[0][i][j][n])
                            # Winner takes the position he wanted (he's the only remaining here)
                            pedRes[0][i][j] = [pedRes[0][i][j][chosenIndex]]
                            pedRes[1][i][j] = [pedRes[1][i][j][chosenIndex]]
                            # Losers are moved back to their original position
                            for n in range(len(remPed)):
                                pedID = remPed[n]
                                I, J = pedList[pedID][0], pedList[pedID][1]
                                pedRes[0][I][J] = [pedID]
                                pedRes[1][I][J] = [[0,0]]
                
                # Move pedestrian and update position, direction and dynamic floor field
                pedPos, pedMoved = np.zeros((iSize,jSize)), np.zeros((iSize,jSize))
                for i in range(iSize):
                    for j in range(jSize):
                        for n in range(len(pedRes[0][i][j])):
                            pedID = pedRes[0][i][j][n]
                            I, J = pedList[pedID][0], pedList[pedID][1]
                            pedList[pedID] = [i,j]
                            pedDir[pedID] = pedRes[1][i][j][n]
                            if i!=I or j!=J:
                                dynamicFF[I][J] = dynamicFF[I][J]+1
                        if len(pedRes[0][i][j])>0:
                            pedPos[i][j] = 1
                            if abs(pedDir[pedID][0])+abs(pedDir[pedID][1])>0:
                                pedMoved[i][j] = 1

                # Store data for current iteration
                simTime = simTime + (self.meshSize / self.vFree)
                if sum(sum(pedPos[0:iSize-3, 0:iSize-3])) > 0 and simTime <= maxTime:
                    if tZero == 0:        tZero = simTime
                    repSimulationPos.append(copy.deepcopy(pedPos[0:iSize-3, 0:iSize-3]))
                    repSimulationTime.append(copy.deepcopy(simTime - tZero))
                    repSimulationHeatmap = repSimulationHeatmap+pedMoved[0:iSize-3, 0:iSize-3] * (self.meshSize / self.vFree)       
            
            #We put repSimulationPos, repSimulationTime and repSimulationHeatmap in a list. It's of size (1, 2+p^2+nt+(p^2)*nt)
            #Notice we would expect all datasets simulated or observed should be stored like this and while computing distance
            #we will use this convention. - This is done as ABCpy expects a dataset as a nparray.
            nTime = math.ceil(self.maxSimulationTime / ((self.meshSize / self.vFree)))
            nCells = np.array(repSimulationHeatmap).shape[0]
            totElements = 2 + nCells**2 + nTime + (nCells**2) * nTime
            # Adding primary information
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

class FFSimulator2(ProbabilisticModel, Continuous):
    """
        Simulation based on the model by Li et al.
        Block-based floor field model for pedestrian’s walking through corner
        Shengnan Li, Xingang Li. Yunchao Qu, Bin Jia
    
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
        self.totalPed = 42                  # number of pedestrians (whole group)
        self.vFree = 1.4                    # free walking speed (m/s)
        self.startSurface = 16              # starting area (m2)

        # Internal settings and options (determine details of the simulation)
        self.motionType = 2             # rule used to determine motion (1 = deterministic approach, 2 = probabilistic approach)
        self.entranceType = 1           # method used to introduce pedestrians (1 = set flow, 2 = use waiting area)
        self.meshSize = 0.4             # mesh size (m)
        self.allowStop = 0              # allow stopping (if 0 stopping not allowed)
        self.maxInflowIter = 100        # iteration limits for inflow (avoid being stuck in a while loop)
        self.maxSimulationTime = 60     # time limit for simulation (s) (in case simulation get stuck due to odd parameter's choice)

        # We expect input of type parameters = [theta1, theta2, n_timestep]
        if not isinstance(parameters, list):
            raise TypeError('Input of FloorField model is of type list')

        if len(parameters) != 7:
            raise RuntimeError('Input list must be of length 7, containing [kS, kD, kW, decay, diffusion, liBufferRatio, liBufferAngle].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)


    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 7:
            raise ValueError('Number of parameters of FloorField model must be 7.')

        # Check whether input is from correct domain
        kS = input_values[0]                # wall floor field weighting factor
        kD = input_values[1]                # static floor field weighting factor
        kW = input_values[2]                # dynamic floor field weighting factor
        decay = input_values[3]             # decay of the dynamic FF
        diffusion = input_values[4]         # diffusion of the dynamic FF
        liBufferRatio = input_values[5]     # length ratio of buffer region
        liBufferAngle = input_values[6]     # deflection angle in the buffer region (deg) [45,90)

        if not isinstance(kS, (float, np.float64, np.float32, np.float16)) or kS < 0:
            print('kS is not of correct type or out of range')
            return False

        if not isinstance(kD, (float, np.float64, np.float32, np.float16)) or kD < 0:
            print('kD is not of correct type or out of range')
            return False
        
        if not isinstance(kW, (float, np.float64, np.float32, np.float16)) or kW < 0:
            print('kW is not of correct type or out of range')
            return False

        if not isinstance(decay, (float, np.float64, np.float32, np.float16)) or decay < 0 or diffusion > 1:
            print('decay is not of correct type or out of range')
            return False

        if not isinstance(diffusion, (float, np.float64, np.float32, np.float16)) or diffusion < 0 or diffusion > 1:
            print('diffusion is not of correct type or out of range')
            return False

        if not isinstance(liBufferRatio, (float, np.float64, np.float32, np.float16)) or liBufferRatio < 0:
            print('liBufferRatio is not of correct type or out of range')
            return False

        if not isinstance(liBufferAngle, (float, np.float64, np.float32, np.float16)) or liBufferAngle < 45 or liBufferAngle > 90:
            print('liBufferAngle is not of correct type or out of range')
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
        kS = input_values[0]                # wall floor field weighting factor
        kD = input_values[1]                # static floor field weighting factor
        kW = input_values[2]                # dynamic floor field weighting factor
        decay = input_values[3]             # decay of the dynamic FF
        diffusion = input_values[4]         # diffusion of the dynamic FF
        liBufferRatio = input_values[5]       # length ratio of buffer region
        liBufferAngle = input_values[6]       # deflection angle in the buffer region (deg) [45,90)
        
        # Do the actual forward simulation
        vector_of_k_samples = self.floorfield(kS, kD, kW, decay, diffusion, liBufferRatio, liBufferAngle, k)
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result


    def floorfield(self, kS, kD, kW, decay, diffusion, liBufferRatio, liBufferAngle, simRep):
        # Generate n_simulated data
        result = []

        ## CREATE STATIC AND WALL FLOOR FIELDS (both depends on parameters chosen)
        # In floorMap: 0 is wall or obstacles       -1 is entrance       2 is a special region not considered in results
        #              1 is walkable area           -2 is exit             (it is used to directly push people to the exit)

        # Read floor map and determine section dimension and corridor width
        if self.entranceType==1:
            floorMap = np.array(list(csv.reader(open("floor_map/Lcorner.csv"), delimiter=","))).astype("float")
            (iSize, jSize) = floorMap.shape
            entranceLocation = [math.inf, 0]
            for i in range(iSize):
                if floorMap[i][jSize - 1] == -1:
                    entranceLocation[0] = min([i, entranceLocation[0]]);
                    entranceLocation[1] = max([i, entranceLocation[1]])
            corrWidth = entranceLocation[1] - entranceLocation[0] + 1
            
        if self.entranceType==2:
            floorMap = np.array(list(csv.reader(open("floor_map/Lcorner_long.csv"), delimiter=","))).astype("float")
            (iSize,jSize) = floorMap.shape
            corrWidth = int(abs(sum(floorMap[:,-1])))
            startArea = [[iSize,0],[jSize,0]]
            for i in range(iSize):
                for j in range(jSize):
                    if floorMap[i,j]==-1:
                        if startArea[0][0]>i:  startArea[0][0] = i
                        if startArea[0][1]<i:  startArea[0][1] = i
                        if startArea[1][0]>j:  startArea[1][0] = j
                        if startArea[1][1]<j:  startArea[1][1] = j

        # Determine location of turning point
        maxDistCenter, iCenter, jCenter = 0, 0, 0
        for i in range(iSize - 1, 0, -1):
            for j in range(jSize - 1, 0, -1):
                dist = (iSize - i) ** 2 + (jSize - j) ** 2
                if dist > maxDistCenter and floorMap[i][j] == 0:        maxDistCenter, iCenter, jCenter = dist, i, j
        
        # Static Floor Field using the block approach by Li et al.           
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
        bufferCells, slope = round(liBufferRatio * corrWidth), math.tan(abs(liBufferAngle - 90) * (math.pi / 180))

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

        minExp, maxExp = math.log(sys.float_info.min*5), math.log(sys.float_info.max/5)

        def compProb(sFFcenter, sFF, dFF, wFF):
            exponent = -kS * (sFF - sFFcenter) + kD * dFF + kW * wFF
            if exponent > minExp and exponent < maxExp:
                return math.exp(-kS * (sFF - sFFcenter) + kD * dFF + kW * wFF)
            elif exponent < minExp:
                return math.exp(minExp)
            elif exponent > maxExp:
                return math.exp(maxExp)

        compProb = np.vectorize(compProb)

        # Main simulation loop  (including repetitions)
        for repi in range(simRep):
            
            # Initialize variables (relative to single execution)
            repSimulationPos, repSimulationTime = [], []
            pedPos, repSimulationHeatmap = np.zeros((iSize, jSize)), np.zeros((iSize-3,iSize-3))
            dynamicFF = np.zeros((iSize, jSize))
            pedList, pedDir = [], []
            numPed, inTime, simTime, tZero = 0, 0, 0, 0
            
            # Place pedestrians in starting area
            if self.entranceType==2:
                while (numPed < self.totalPed):
                    m, n = rand.randint(startArea[0][0],startArea[0][1]), rand.randint(startArea[1][0],startArea[1][1])
                    if pedPos[m,n]==0:
                        pedPos[m][n] = 1
                        numPed = numPed + 1
                        pedList.append([m,n])
                        pedDir.append([0,-1])

            # Internal simulation loop (single run)
            maxTime = self.maxSimulationTime - (self.meshSize/self.vFree)
            while ((numPed < self.totalPed) or (len(pedList) > 0)) and (simTime-inTime < maxTime):
                
                # Inject pedestrians if required (entrance by inflow value)
                if self.entranceType==1:
                    if numPed < self.totalPed:
                        warningCounter = 0
                        while (numPed < self.vFree * corrWidth * self.meshSize * (self.totalPed / self.startSurface) * (simTime - inTime)):
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

                 # Define starting time for starting area setup
                if self.entranceType==2:
                    if sum(pedPos[:,iSize])>1 and inTime==0:
                        inTime = simTime                
                
                # Remove pedestrians if required
                remPed = []
                for pedID in range(len(pedList)):
                    if abs(floorMap[pedList[pedID][0]][pedList[pedID][1]]) == 2:      remPed.append(pedID)
                remPed.sort(reverse=True)
                for pedID in range(len(remPed)):
                    del pedList[remPed[pedID]]
                    del pedDir[remPed[pedID]]
                    
                # Diffuse and decay dynamic floor field
                for i in range(iSize):
                    for j in range(jSize):
                        if floorMap[i][j] != 0:
                            # Decay
                            localDynamicFF = int(dynamicFF[i][j])
                            if dynamicFF[i][j] > 0:
                                for n in range(0,localDynamicFF):
                                    if rand.uniform(0,1) < decay:
                                        dynamicFF[i][j] = dynamicFF[i][j] - 1
                            # Diffusion
                            if dynamicFF[i][j] > 0 and rand.uniform(0,1) < diffusion:
                                dynamicFF[i][j] = dynamicFF[i][j] - 1
                                index = rand.randint(0,1)
                                if index==0:
                                    dI, dJ = int(np.sign(rand.uniform(-1,1))), 0
                                elif index==1:
                                    dI, dJ = 0, int(np.sign(rand.uniform(-1,1)))
                                if i+dI < iSize and j+dJ < jSize and floorMap[i+dI][j+dJ]!=0:
                                    dynamicFF[i+dI][j+dJ] = dynamicFF[i+dI][j+dJ] + 1    

                # Compute transition probabilities and reserve potential positions
                pedRes = [[[[] for j in range(jSize)] for i in range(iSize)] for n in range(2)]
                for pedID in range(len(pedList)):
                    
                    # Extract local floor fields and prepare data for transition probability calculation
                    I, dI, J, dJ = pedList[pedID][0], pedDir[pedID][0], pedList[pedID][1], pedDir[pedID][1]
                    localMap = [None]*3
                    localMap[0] = compProb(staticFF[I,J], staticFF[I-1:I+2,J-1:J+2], dynamicFF[I-1:I+2,J-1:J+2], wallFF[I-1:I+2,J-1:J+2])
                    localMap[1], localMap[2] = pedPos[I-1:I+2,J-1:J+2], floorMap[I-1:I+2,J-1:J+2]
                    if (dI!=0 or dJ!=0) and I-dI>=0 and I-dI<iSize and J-dJ>=0 and J-dJ<jSize:
                        localValue = compProb(staticFF[I,J], staticFF[I-dI,J-dJ], dynamicFF[I-dI,J-dJ]-1, wallFF[I-dI,J-dJ])
                        localMap[0][1-dI][1-dJ] = localValue.item(0)
                    for n in range(len(localMap)):
                        if J==jSize-1:      localMap[n] = np.c_[localMap[n], np.zeros(3)]
                        elif J==0:          localMap[n] = np.c_[np.zeros(3), localMap[n]]
                        elif I==iSize-1:    localMap[n] = np.r_[localMap[n], [np.zeros(3)]]
                        elif I==0:          localMap[n] = np.r_[[np.zeros(3)], localMap[n]]
                        
                    # Exclude corners (Von Neumann neighborhood), back step, stop (if not allowed) and inacessible locations
                    localMap[0][0][0], localMap[0][2][0], localMap[0][0][2], localMap[0][2][2] = 0, 0, 0, 0
                    
                    localMap[0][1-dI][1-dJ] = 0                             # remove back stepping
                    if self.allowStop==0:        localMap[0][1][1] = 0      # no stopping allowed (depending on options)
                    for i in range(localMap[0].shape[0]):
                        for j in range(localMap[0].shape[1]):
                            if not(i==1 and j==1) and localMap[1][i][j]==1:     localMap[0][i][j] = 0       # exclude occupied positions
                            if localMap[2][i][j]==0:                            localMap[0][i][j] = 0       # exclude inacessible locations   
                    totalValue = sum(sum(localMap[0]))
                    
                    # Determine cell to reserve for next step
                    if self.motionType==1:
                        # Determine maximum values and relative index (also for tied cases)
                        if totalValue>0:    transProb = localMap[0]/totalValue          # compute transition probability
                        maxValue, maxIndex = transProb.max(), []
                        if totalValue>0:
                            for i in range(transProb.shape[0]):
                                for j in range(transProb.shape[1]):
                                    if transProb[i][j]==maxValue:
                                        maxIndex.append([i,j])
                            if len(maxIndex)==1:        m, n = maxIndex[0]      # only one maximum
                            else:                       m, n = maxIndex[rand.randint(0,len(maxIndex)-1)]
                        else:   m, n = 1, 1         # special case where transition probability is 0 everywhere
                    elif self.motionType==2:
                        # Determine next cell where do move stochastically based on previous results
                        p0, p1, p2, p3, p4 = localMap[0][1][1], localMap[0][1][0], localMap[0][1][2], localMap[0][0][1], localMap[0][2][1]
                        if totalValue>0:
                            randValue = rand.uniform(0,1)
                            if randValue >= (p0/totalValue) and randValue < ((p0+p1)/totalValue):
                                m, n = 1, 0
                            elif randValue >= ((p0+p1)/totalValue) and randValue < ((p0+p1+p2)/totalValue):
                                m, n = 1, 2
                            elif randValue >= ((p0+p1+p2)/totalValue) and randValue < ((p0+p1+p2+p3)/totalValue):
                                m, n = 0, 1
                            elif randValue >= ((p0+p1+p2+p3)/totalValue) and randValue < ((p0+p1+p2+p3+p4)/totalValue):
                                m, n = 2, 1
                            else:
                                m, n = 1, 1
                        else:
                            m, n = 1, 1
                        
                    # Reserve position
                    if I+(m-1)>=0 and I+(m-1)<iSize and J+(n-1)>=0 and J+(n-1)<jSize:
                        pedRes[0][I+(m-1)][J+(n-1)].append(pedID)
                        if floorMap[I+(m-1)][J+(n-1)]==2:   pedRes[1][I+(m-1)][J+(n-1)].append([dI,dJ])
                        else:                               pedRes[1][I+(m-1)][J+(n-1)].append([m-1,n-1])
                    
                # Resolve conflicts (in areas with floorMap==2, conflics are possible to avoid congestion near exit)
                for i in range(iSize):
                    for j in range(jSize-2):
                        nReserved = len(pedRes[0][i][j])
                        if nReserved>1:
                            # Determine winner candidate
                            chosenIndex, remPed = rand.randint(0,nReserved-1), []
                            # Losers are listed to move them back later
                            for n in range(nReserved):
                                if n!=chosenIndex:      remPed.append(pedRes[0][i][j][n])
                            # Winner takes the position he wanted (he's the only remaining here)
                            pedRes[0][i][j] = [pedRes[0][i][j][chosenIndex]]
                            pedRes[1][i][j] = [pedRes[1][i][j][chosenIndex]]
                            # Losers are moved back to their original position
                            for n in range(len(remPed)):
                                pedID = remPed[n]
                                I, J = pedList[pedID][0], pedList[pedID][1]
                                pedRes[0][I][J] = [pedID]
                                pedRes[1][I][J] = [[0,0]]
                
                # Move pedestrian and update position, direction and dynamic floor field
                pedPos, pedMoved = np.zeros((iSize,jSize)), np.zeros((iSize,jSize))
                for i in range(iSize):
                    for j in range(jSize):
                        for n in range(len(pedRes[0][i][j])):
                            pedID = pedRes[0][i][j][n]
                            I, J = pedList[pedID][0], pedList[pedID][1]
                            pedList[pedID] = [i,j]
                            pedDir[pedID] = pedRes[1][i][j][n]
                            if i!=I or j!=J:
                                dynamicFF[I][J] = dynamicFF[I][J]+1
                        if len(pedRes[0][i][j])>0:
                            pedPos[i][j] = 1
                            if abs(pedDir[pedID][0])+abs(pedDir[pedID][1])>0:
                                pedMoved[i][j] = 1

                # Store data for current iteration
                simTime = simTime + (self.meshSize / self.vFree)
                if sum(sum(pedPos[0:iSize-3, 0:iSize-3])) > 0 and simTime <= maxTime:
                    if tZero == 0:        tZero = simTime
                    repSimulationPos.append(copy.deepcopy(pedPos[0:iSize-3, 0:iSize-3]))
                    repSimulationTime.append(copy.deepcopy(simTime - tZero))
                    repSimulationHeatmap = repSimulationHeatmap+pedMoved[0:iSize-3, 0:iSize-3] * (self.meshSize / self.vFree)       
            
            #We put repSimulationPos, repSimulationTime and repSimulationHeatmap in a list. It's of size (1, 2+p^2+nt+(p^2)*nt)
            #Notice we would expect all datasets simulated or observed should be stored like this and while computing distance
            #we will use this convention. - This is done as ABCpy expects a dataset as a nparray.
            nTime = math.ceil(self.maxSimulationTime / ((self.meshSize / self.vFree)))
            nCells = np.array(repSimulationHeatmap).shape[0]
            totElements = 2 + nCells**2 + nTime + (nCells**2) * nTime
            # Adding primary information
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

class FFSimulator3(ProbabilisticModel, Continuous):
    """
        Simulation based on the model by Dias and Lovreglio
        Calibrating cellular automaton models for pedestrians walking through corners
        Charitha Dias, Ruggiero Lovreglio
    
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
        self.totalPed = 42                  # number of pedestrians (whole group)
        self.vFree = 1.4                    # free walking speed (m/s)
        self.startSurface = 16              # starting area (m2)

        # Internal settings and options (determine details of the simulation)
        self.motionType = 2             # rule used to determine motion (1 = deterministic approach, 2 = probabilistic approach)
        self.entranceType = 1           # method used to introduce pedestrians (1 = set flow, 2 = use waiting area)
        self.meshSize = 0.4             # mesh size (m)
        self.allowStop = 0              # allow stopping (if 0 stopping not allowed)
        self.maxInflowIter = 100        # iteration limits for inflow (avoid being stuck in a while loop)
        self.maxSimulationTime = 60     # time limit for simulation (s) (in case simulation get stuck due to odd parameter's choice)

        # We expect input of type parameters = [theta1, theta2, n_timestep]
        if not isinstance(parameters, list):
            raise TypeError('Input of FloorField model is of type list')

        if len(parameters) != 9:
            raise RuntimeError('Input list must be of length 9, containing [kS, kD, kW, decay, diffusion, diasA, diasB, diasP1, diasP2].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)


    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 9:
            raise ValueError('Number of parameters of FloorField model must be 9.')

        # Check whether input is from correct domain
        kS = input_values[0]                # wall floor field weighting factor
        kD = input_values[1]                # static floor field weighting factor
        kW = input_values[2]                # dynamic floor field weighting factor
        decay = input_values[3]             # decay of the dynamic FF
        diffusion = input_values[4]         # diffusion of the dynamic FF
        diasA = input_values[5]             # turning point distance a in Dias model
        diasB = input_values[6]             # turning point distance b in Dias model
        diasP1 = input_values[7]            # distribution p1
        diasP2 = input_values[8]            # distribution p2

        if not isinstance(kS, (float, np.float64, np.float32, np.float16)) or kS < 0:
            print('kS is not of correct type or out of range')
            return False

        if not isinstance(kD, (float, np.float64, np.float32, np.float16)) or kD < 0:
            print('kD is not of correct type or out of range')
            return False
        
        if not isinstance(kW, (float, np.float64, np.float32, np.float16)) or kW < 0:
            print('kW is not of correct type or out of range')
            return False

        if not isinstance(decay, (float, np.float64, np.float32, np.float16)) or decay < 0 or diffusion > 1:
            print('decay is not of correct type or out of range')
            return False

        if not isinstance(diffusion, (float, np.float64, np.float32, np.float16)) or diffusion < 0 or diffusion > 1:
            print('diffusion is not of correct type or out of range')
            return False

        if not isinstance(diasA, (float, np.float64, np.float32, np.float16)) or diasA < 0:
            print('diasA is not of correct type or out of range')
            return False

        if not isinstance(diasB, (float, np.float64, np.float32, np.float16)) or diasB < 0:
            print('diasB is not of correct type or out of range')
            return False
        
        if not isinstance(diasP1, (float, np.float64, np.float32, np.float16)) or diasP1 < 0 or diasP1 > 1:
            print('diasP1 is not of correct type or out of range')
            return False

        if not isinstance(diasP2, (float, np.float64, np.float32, np.float16)) or diasP2 < 0 or diasP2 > 1:
            print('diasP2 is not of correct type or out of range')
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
        kS = input_values[0]                # wall floor field weighting factor
        kD = input_values[1]                # static floor field weighting factor
        kW = input_values[2]                # dynamic floor field weighting factor
        decay = input_values[3]             # decay of the dynamic FF
        diffusion = input_values[4]         # diffusion of the dynamic FF
        diasA = input_values[5]             # turning point distance a in Dias model
        diasB = input_values[6]             # turning point distance b in Dias model
        diasP1 = input_values[7]            # distribution p1
        diasP2 = input_values[8]            # distribution p2
        
        # Do the actual forward simulation
        vector_of_k_samples = self.floorfield(kS, kD, kW, decay, diffusion, diasA, diasB, diasP1, diasP2, k)
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result


    def floorfield(self, kS, kD, kW, decay, diffusion, diasA, diasB, diasP1, diasP2, simRep):
        # Generate n_simulated data
        result = []

        ## CREATE STATIC AND WALL FLOOR FIELDS (both depends on parameters chosen)
        # In floorMap: 0 is wall or obstacles       -1 is entrance       2 is a special region not considered in results
        #              1 is walkable area           -2 is exit             (it is used to directly push people to the exit)

        # Read floor map and determine section dimension and corridor width
        if self.entranceType==1:
            floorMap = np.array(list(csv.reader(open("floor_map/Lcorner.csv"), delimiter=","))).astype("float")
            (iSize, jSize) = floorMap.shape
            entranceLocation = [math.inf, 0]
            for i in range(iSize):
                if floorMap[i][jSize - 1] == -1:
                    entranceLocation[0] = min([i, entranceLocation[0]]);
                    entranceLocation[1] = max([i, entranceLocation[1]])
            corrWidth = entranceLocation[1] - entranceLocation[0] + 1
            
        if self.entranceType==2:
            floorMap = np.array(list(csv.reader(open("floor_map/Lcorner_long.csv"), delimiter=","))).astype("float")
            (iSize,jSize) = floorMap.shape
            corrWidth = int(abs(sum(floorMap[:,-1])))
            startArea = [[iSize,0],[jSize,0]]
            for i in range(iSize):
                for j in range(jSize):
                    if floorMap[i,j]==-1:
                        if startArea[0][0]>i:  startArea[0][0] = i
                        if startArea[0][1]<i:  startArea[0][1] = i
                        if startArea[1][0]>j:  startArea[1][0] = j
                        if startArea[1][1]<j:  startArea[1][1] = j

        # Determine location of turning point
        maxDistCenter, iCenter, jCenter = 0, 0, 0
        for i in range(iSize - 1, 0, -1):
            for j in range(jSize - 1, 0, -1):
                dist = (iSize - i) ** 2 + (jSize - j) ** 2
                if dist > maxDistCenter and floorMap[i][j] == 0:        maxDistCenter, iCenter, jCenter = dist, i, j
        
        # Define turning center and function for turning angle
        iTurning, jTurning = iCenter + math.floor(diasB / self.meshSize) - 1, jCenter + math.floor(diasA / self.meshSize) - 1
        xTurning, yTurning = iCenter + (diasB / self.meshSize) - 1, jCenter + (diasA / self.meshSize) - 1
        def alpha(theta):
            return beta.cdf(theta / (math.pi / 2), diasP1, diasP2) * (math.pi / 2)
        
        # Create floor map for last straight section
        staticFF = np.zeros((iSize,jSize))
        for i in range(iSize-1,iTurning,-1):
            if i>=1:        staticFF[i,1:corrWidth+1] = iSize-i-1   
        
        # Create static FF for turning part
        alphaMatrix, maxStaticFF = np.zeros((iSize,jSize)), np.amax(staticFF)
        for i in range(0, iTurning + 1):
            for j in range(0, jTurning + 1):
                if floorMap[i,j]==1:
                    iDist, jDist = xTurning - i + 0.5, yTurning - j + 0.5
                    if iDist==0:    alphaMatrix[i,j] = math.degrees(alpha(math.pi / 2))
                    else:           alphaMatrix[i,j] = math.degrees(alpha(math.atan(jDist/iDist)))
                    alphaMatrix[i,j] = maxStaticFF + (1 - (alphaMatrix[i,j] / 90))
        staticFF = alphaMatrix + staticFF 
        
        # Create floor map for initial approach (initial straight section)
        startValue = math.ceil(np.amax(staticFF)) + 1
        for j in range(jTurning + 1, jSize):
            staticFF[1:corrWidth+1,j] = startValue + (j - jTurning - 1)
    
        # Normalize to have only integer values
        allData = staticFF.flatten()
        setValues = sorted(set(allData))
        for i in range(iSize):
            for j in range(jSize):
                staticFF[i,j] = setValues.index(staticFF[i,j])

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

        minExp, maxExp = math.log(sys.float_info.min*5), math.log(sys.float_info.max/5)

        def compProb(sFFcenter, sFF, dFF, wFF):
            exponent = -kS * (sFF - sFFcenter) + kD * dFF + kW * wFF
            if exponent > minExp and exponent < maxExp:
                return math.exp(-kS * (sFF - sFFcenter) + kD * dFF + kW * wFF)
            elif exponent < minExp:
                return math.exp(minExp)
            elif exponent > maxExp:
                return math.exp(maxExp)

        compProb = np.vectorize(compProb)

        # Main simulation loop  (including repetitions)
        for repi in range(simRep):
            
            # Initialize variables (relative to single execution)
            repSimulationPos, repSimulationTime = [], []
            pedPos, repSimulationHeatmap = np.zeros((iSize, jSize)), np.zeros((iSize-3,iSize-3))
            dynamicFF = np.zeros((iSize, jSize))
            pedList, pedDir = [], []
            numPed, inTime, simTime, tZero = 0, 0, 0, 0
            
            # Place pedestrians in starting area
            if self.entranceType==2:
                while (numPed < self.totalPed):
                    m, n = rand.randint(startArea[0][0],startArea[0][1]), rand.randint(startArea[1][0],startArea[1][1])
                    if pedPos[m,n]==0:
                        pedPos[m][n] = 1
                        numPed = numPed + 1
                        pedList.append([m,n])
                        pedDir.append([0,-1])

            # Internal simulation loop (single run)
            maxTime = self.maxSimulationTime - (self.meshSize/self.vFree)
            while ((numPed < self.totalPed) or (len(pedList) > 0)) and (simTime-inTime < maxTime):
                
                # Inject pedestrians if required (entrance by inflow value)
                if self.entranceType==1:
                    if numPed < self.totalPed:
                        warningCounter = 0
                        while (numPed < self.vFree * corrWidth * self.meshSize * (self.totalPed / self.startSurface) * (simTime - inTime)):
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

                 # Define starting time for starting area setup
                if self.entranceType==2:
                    if sum(pedPos[:,iSize])>1 and inTime==0:
                        inTime = simTime                
                
                # Remove pedestrians if required
                remPed = []
                for pedID in range(len(pedList)):
                    if abs(floorMap[pedList[pedID][0]][pedList[pedID][1]]) == 2:      remPed.append(pedID)
                remPed.sort(reverse=True)
                for pedID in range(len(remPed)):
                    del pedList[remPed[pedID]]
                    del pedDir[remPed[pedID]]
                    
                # Diffuse and decay dynamic floor field
                for i in range(iSize):
                    for j in range(jSize):
                        if floorMap[i][j] != 0:
                            # Decay
                            localDynamicFF = int(dynamicFF[i][j])
                            if dynamicFF[i][j] > 0:
                                for n in range(0,localDynamicFF):
                                    if rand.uniform(0,1) < decay:
                                        dynamicFF[i][j] = dynamicFF[i][j] - 1
                            # Diffusion
                            if dynamicFF[i][j] > 0 and rand.uniform(0,1) < diffusion:
                                dynamicFF[i][j] = dynamicFF[i][j] - 1
                                index = rand.randint(0,1)
                                if index==0:
                                    dI, dJ = int(np.sign(rand.uniform(-1,1))), 0
                                elif index==1:
                                    dI, dJ = 0, int(np.sign(rand.uniform(-1,1)))
                                if i+dI < iSize and j+dJ < jSize and floorMap[i+dI][j+dJ]!=0:
                                    dynamicFF[i+dI][j+dJ] = dynamicFF[i+dI][j+dJ] + 1    

                # Compute transition probabilities and reserve potential positions
                pedRes = [[[[] for j in range(jSize)] for i in range(iSize)] for n in range(2)]
                for pedID in range(len(pedList)):
                    
                    # Extract local floor fields and prepare data for transition probability calculation
                    I, dI, J, dJ = pedList[pedID][0], pedDir[pedID][0], pedList[pedID][1], pedDir[pedID][1]
                    localMap = [None]*3
                    localMap[0] = compProb(staticFF[I,J], staticFF[I-1:I+2,J-1:J+2], dynamicFF[I-1:I+2,J-1:J+2], wallFF[I-1:I+2,J-1:J+2])
                    localMap[1], localMap[2] = pedPos[I-1:I+2,J-1:J+2], floorMap[I-1:I+2,J-1:J+2]
                    if (dI!=0 or dJ!=0) and I-dI>=0 and I-dI<iSize and J-dJ>=0 and J-dJ<jSize:
                        localValue = compProb(staticFF[I,J], staticFF[I-dI,J-dJ], dynamicFF[I-dI,J-dJ]-1, wallFF[I-dI,J-dJ])
                        localMap[0][1-dI][1-dJ] = localValue.item(0)
                    for n in range(len(localMap)):
                        if J==jSize-1:      localMap[n] = np.c_[localMap[n], np.zeros(3)]
                        elif J==0:          localMap[n] = np.c_[np.zeros(3), localMap[n]]
                        elif I==iSize-1:    localMap[n] = np.r_[localMap[n], [np.zeros(3)]]
                        elif I==0:          localMap[n] = np.r_[[np.zeros(3)], localMap[n]]
                        
                    # Exclude corners (Von Neumann neighborhood), back step, stop (if not allowed) and inacessible locations
                    localMap[0][0][0], localMap[0][2][0], localMap[0][0][2], localMap[0][2][2] = 0, 0, 0, 0
                    
                    localMap[0][1-dI][1-dJ] = 0                             # remove back stepping
                    if self.allowStop==0:        localMap[0][1][1] = 0      # no stopping allowed (depending on options)
                    for i in range(localMap[0].shape[0]):
                        for j in range(localMap[0].shape[1]):
                            if not(i==1 and j==1) and localMap[1][i][j]==1:     localMap[0][i][j] = 0       # exclude occupied positions
                            if localMap[2][i][j]==0:                            localMap[0][i][j] = 0       # exclude inacessible locations   
                    totalValue = sum(sum(localMap[0]))
                    
                    # Determine cell to reserve for next step
                    if self.motionType==1:
                        # Determine maximum values and relative index (also for tied cases)
                        if totalValue>0:    transProb = localMap[0]/totalValue          # compute transition probability
                        maxValue, maxIndex = transProb.max(), []
                        if totalValue>0:
                            for i in range(transProb.shape[0]):
                                for j in range(transProb.shape[1]):
                                    if transProb[i][j]==maxValue:
                                        maxIndex.append([i,j])
                            if len(maxIndex)==1:        m, n = maxIndex[0]      # only one maximum
                            else:                       m, n = maxIndex[rand.randint(0,len(maxIndex)-1)]
                        else:   m, n = 1, 1         # special case where transition probability is 0 everywhere
                    elif self.motionType==2:
                        # Determine next cell where do move stochastically based on previous results
                        p0, p1, p2, p3, p4 = localMap[0][1][1], localMap[0][1][0], localMap[0][1][2], localMap[0][0][1], localMap[0][2][1]
                        if totalValue>0:
                            randValue = rand.uniform(0,1)
                            if randValue >= (p0/totalValue) and randValue < ((p0+p1)/totalValue):
                                m, n = 1, 0
                            elif randValue >= ((p0+p1)/totalValue) and randValue < ((p0+p1+p2)/totalValue):
                                m, n = 1, 2
                            elif randValue >= ((p0+p1+p2)/totalValue) and randValue < ((p0+p1+p2+p3)/totalValue):
                                m, n = 0, 1
                            elif randValue >= ((p0+p1+p2+p3)/totalValue) and randValue < ((p0+p1+p2+p3+p4)/totalValue):
                                m, n = 2, 1
                            else:
                                m, n = 1, 1
                        else:
                            m, n = 1, 1
                        
                    # Reserve position
                    if I+(m-1)>=0 and I+(m-1)<iSize and J+(n-1)>=0 and J+(n-1)<jSize:
                        pedRes[0][I+(m-1)][J+(n-1)].append(pedID)
                        if floorMap[I+(m-1)][J+(n-1)]==2:   pedRes[1][I+(m-1)][J+(n-1)].append([dI,dJ])
                        else:                               pedRes[1][I+(m-1)][J+(n-1)].append([m-1,n-1])
                    
                # Resolve conflicts (in areas with floorMap==2, conflics are possible to avoid congestion near exit)
                for i in range(iSize):
                    for j in range(jSize-2):
                        nReserved = len(pedRes[0][i][j])
                        if nReserved>1:
                            # Determine winner candidate
                            chosenIndex, remPed = rand.randint(0,nReserved-1), []
                            # Losers are listed to move them back later
                            for n in range(nReserved):
                                if n!=chosenIndex:      remPed.append(pedRes[0][i][j][n])
                            # Winner takes the position he wanted (he's the only remaining here)
                            pedRes[0][i][j] = [pedRes[0][i][j][chosenIndex]]
                            pedRes[1][i][j] = [pedRes[1][i][j][chosenIndex]]
                            # Losers are moved back to their original position
                            for n in range(len(remPed)):
                                pedID = remPed[n]
                                I, J = pedList[pedID][0], pedList[pedID][1]
                                pedRes[0][I][J] = [pedID]
                                pedRes[1][I][J] = [[0,0]]
                
                # Move pedestrian and update position, direction and dynamic floor field
                pedPos, pedMoved = np.zeros((iSize,jSize)), np.zeros((iSize,jSize))
                for i in range(iSize):
                    for j in range(jSize):
                        for n in range(len(pedRes[0][i][j])):
                            pedID = pedRes[0][i][j][n]
                            I, J = pedList[pedID][0], pedList[pedID][1]
                            pedList[pedID] = [i,j]
                            pedDir[pedID] = pedRes[1][i][j][n]
                            if i!=I or j!=J:
                                dynamicFF[I][J] = dynamicFF[I][J]+1
                        if len(pedRes[0][i][j])>0:
                            pedPos[i][j] = 1
                            if abs(pedDir[pedID][0])+abs(pedDir[pedID][1])>0:
                                pedMoved[i][j] = 1

                # Store data for current iteration
                simTime = simTime + (self.meshSize / self.vFree)
                if sum(sum(pedPos[0:iSize-3, 0:iSize-3])) > 0 and simTime <= maxTime:
                    if tZero == 0:        tZero = simTime
                    repSimulationPos.append(copy.deepcopy(pedPos[0:iSize-3, 0:iSize-3]))
                    repSimulationTime.append(copy.deepcopy(simTime - tZero))
                    repSimulationHeatmap = repSimulationHeatmap+pedMoved[0:iSize-3, 0:iSize-3] * (self.meshSize / self.vFree)       
            
            #We put repSimulationPos, repSimulationTime and repSimulationHeatmap in a list. It's of size (1, 2+p^2+nt+(p^2)*nt)
            #Notice we would expect all datasets simulated or observed should be stored like this and while computing distance
            #we will use this convention. - This is done as ABCpy expects a dataset as a nparray.
            nTime = math.ceil(self.maxSimulationTime / ((self.meshSize / self.vFree)))
            nCells = np.array(repSimulationHeatmap).shape[0]
            totElements = 2 + nCells**2 + nTime + (nCells**2) * nTime
            # Adding primary information
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

