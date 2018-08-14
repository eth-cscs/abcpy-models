import csv
import math
import copy
import numpy as np
from os import listdir
from os.path import isfile, join

def ReadObsdata():
    # Constants (related to experimental or simulation conditions, known a priori)
    maxPed = 16  # number of pedestrians (whole group)
    meshSize = 0.4  # mesh size (m)
    vFree = 1.4  # free walking speed (m/s)
    direction = 'L'  # turning direction
    fps = 29.97  # frame rate in video recordings
    maxTime = 30 # maximum possible time for passing (s)
    corrLength = 2.50  # length of the corridor considered (m)
    repNumber = -1  # number of repetition (-1 loop all)

    ## READ EXPERIMENTAL DATA
    # Read file list and prepare relevant strings
    onlyFiles = [fileList for fileList in listdir('observed_data') if isfile(join('observed_data', fileList))]
    fileClass = '090_UNI_' + direction + format(maxPed, '02d')

    result = []

    # Loop over all data in folder
    for fileName in onlyFiles:
        if fileName.find(fileClass) != -1 and (repNumber == -1 or int(fileName[-5]) == repNumber):

            # Read raw data and create position matrix
            iSize, jSize = math.ceil(corrLength/meshSize)*2-1, math.ceil(corrLength/meshSize)*2-1
            dataRaw = np.array(list(csv.reader(open('observed_data/' + fileName), delimiter=','))).astype('float')
            dataAll = np.zeros((len(dataRaw),7))
            dataAll[:,:-2] = dataRaw
            tMin, tMax = np.min(dataAll[:,2]), np.max(dataAll[:,2])
            timeSteps = math.ceil((tMax - tMin) * fps) + 1
            positionFull = [[] for i in range(timeSteps)]
            for i in range(len(dataAll)):
                localCopy = [dataAll[i,3]*0.01,dataAll[i,4]*0.01]
                positionFull[int(round((dataAll[i,2]-tMin)*fps))].append(localCopy)
                I = math.floor((dataAll[i,3]*0.01+corrLength)/meshSize)
                J = iSize-math.floor((dataAll[i,4]*0.01+corrLength)/meshSize)-1
                dataAll[i,5], dataAll[i,6] = I, J

            # Determine time of first pedestrian in and last out
            tFirst, tLast, repObservationTime = 0, 0, []
            for n in range(timeSteps):
                for i in range(len(positionFull[n])):
                    if abs(positionFull[n][i][0])<corrLength and abs(positionFull[n][i][1])<corrLength:
                        tLast = n
                        if tFirst==0:
                            tFirst = n
            repObservationTime.append(copy.deepcopy((tLast-tFirst)*(1/fps)))
            
            # Delete frames to correspond to simulation conditions
            obsTime = []
            for n in range(tFirst,tLast):
                obsTime.append((n-tFirst)*(1/fps))
            simSteps = math.ceil(max(obsTime)/(meshSize/vFree))
            positionAll = []
            for n in range(simSteps):
                diff, index = 1E+10, 0
                for i in range(len(obsTime)):
                    if abs(n*(meshSize/vFree)-obsTime[i])<diff:
                        diff = abs(n*(meshSize/vFree)-obsTime[i])
                        index = i
                positionAll.append(positionFull[index])
               
            # Create matrix containing positions (including repetitions)
            repObservationPos = np.zeros((iSize,jSize,tLast-tFirst))
            for n in range(len(positionAll)):
                for i in range(len(positionAll[n])):
                    if abs(positionAll[n][i][0])<corrLength and abs(positionAll[n][i][1])<corrLength:
                        I = math.floor((positionAll[n][i][1]+corrLength)/meshSize)
                        J = iSize-math.floor((positionAll[n][i][0]+corrLength)/meshSize)-1
                        repObservationPos[I][J][n-tFirst] = repObservationPos[I][J][n-tFirst]+1   
                        
            # Remove reboundant positions and create related position matrix
            dataUnique = dataAll
            for i in range(len(dataUnique)-1,0,-1):
                if dataUnique[i,5]==dataUnique[i-1,5] and dataUnique[i,6]==dataUnique[i-1,6]:
                    dataUnique = np.delete(dataUnique,i,0)
            positionUnique = [[] for i in range(timeSteps)]
            for i in range(len(dataUnique)):
                localCopy = [dataUnique[i,3]*0.01,dataUnique[i,4]*0.01]
                positionUnique[int(round((dataUnique[i,2]-tMin)*fps))].append(localCopy)
                        
            # Create matrix containing unique positions (exclude stopping cases)
            pedMoved = np.zeros((iSize,jSize,tLast-tFirst))
            for n in range(tFirst,tLast):
                for i in range(len(positionUnique[n])):
                    if abs(positionUnique[n][i][0])<corrLength and abs(positionUnique[n][i][1])<corrLength:
                        I = math.floor((positionUnique[n][i][1]+corrLength)/meshSize)
                        J = iSize-math.floor((positionUnique[n][i][0]+corrLength)/meshSize)-1
                        pedMoved[I][J][n-tFirst] = pedMoved[I][J][n-tFirst]+1     
            
            # Create heatmap
            repObservationHeatmap = np.zeros((iSize-1,jSize-1))
            for n in range(pedMoved.shape[2]):
                repObservationHeatmap = repObservationHeatmap+pedMoved[1:,1:,n]*(meshSize/vFree)
            
            # We put repObservationPos, repObservationTime and repObservationHeatmap in a list. It's of size (1, 2+p^2+(p^2)*nt)
            # Notice we would expect all datasets simulated or observed should be stored like this and while computing distance
            # we will use this convention. - This is done as ABCpy expects a dataset as a nparray.
            nTime = math.ceil(maxTime * fps)
            nCells = np.array(repObservationHeatmap).shape[0]
            totElements = 2 + nCells**2 + nTime + (nCells**2) * nTime
            # Addiing primary information
            repObservation = list([nTime, nCells])
            # Adding the ObservationHeatmap
            repObservation += list(np.array(repObservationHeatmap).reshape(-1,))
            # Adding the ObservationTime
            repObservation += list(repObservationTime)
            for inadd in range(len(repObservation),2 + nCells**2 + nTime):
                repObservation += list(np.array([0]))
            # Adding the ObservationPos
            for indadd in range(len(repObservationPos)):
                repObservation += list(np.array(repObservationPos[indadd]).reshape(-1,))
            # Fill remaining elements
            for indadd in range(len(repObservation),totElements):
                repObservation += list(np.array([0]))

            # Append the result one repetition to the final experimental result
            result.append(copy.deepcopy(repObservation))

    return result
