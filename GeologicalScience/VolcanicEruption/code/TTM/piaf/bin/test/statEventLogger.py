#!/usr/bin/python3

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

def firstline(filename):
    with open(filename, 'r') as f:
        for line in f:
            if(line.split(',')[0].strip() != "time"):
                return line

def stats(filename):
    initTime = int(firstline(filename).split(',')[0].strip())
    lastTime = initTime
    times = []
    accumulatedTime = 0
    with open(filename) as f:
        for line in f:
            if len(line.split(',')) > 0 and line.split(',')[0].strip() != 'time' and int(line.split(',')[0].strip()) != initTime:
                newTime = int(line.split(',')[0].strip())
                accumulatedTime += newTime - lastTime
                times.append(newTime - lastTime)
                lastTime = newTime
    totalTime = lastTime - initTime
    return totalTime - accumulatedTime, times

def stats2(filename):
    with open(filename) as f:
        lines = f.readlines()
    return int(lines[1].split(',')[0]), int(lines[2].split(',')[0]), int(lines[3].split(',')[0])

def trim(data, keep):
    remove = int(len(data)*((1.0-keep)/2.0))
    data.sort()
    return sorted(data, reverse=True)[remove:-remove]

filename = sys.argv[1]
filename2 = sys.argv[2]
timeVariation, times = stats(filename)
keep = 0.97

print("logger statistics of file " + filename + " and " + filename2)
print("difference between total time and cumulated time : " + str(timeVariation) + " [ns]")
print("stats keeping " + str(keep*100) + "% of data :")

print(len(times))
times = trim(times,keep)
print(len(times))

print("mean : " + str(np.mean(times)) + "[ns] == " + str(np.mean(times)/1e6) + "[ms]")
print("std : "+ str(np.std(times)) + "[ns] == " + str(np.std(times)/1e6) + "[ms]")
print("4std/mean : " + str(4*np.std(times)/np.mean(times)))
print("max : "+str(np.max(times)))
print("min : "+str(np.min(times)))
print("mean : "+str(np.mean(times)))
print("(max-min)/mean : " + str((np.max(times)-np.min(times))/np.mean(times)))

times2 = stats2(filename2)
timelog = times2[1]-times2[0]
timenolog = times2[2]-times2[1]

print("execution time without log : " + str(timenolog) + "[ns] == " + str(timenolog/1e6) + "[ms]")
print("execution time with log : " + str(timelog) + "[ns] == " + str(timelog/1e6) + "[ms]")
print("difference between time with and without log : " + str(timelog-timenolog) + "[ns] == " + str((timelog-timenolog)/1e6) + "[ms]")
print("(max-min)/mean : " + str(((np.max([timelog, timenolog])-np.min([timelog, timenolog]))/np.mean([timelog, timenolog]))))

plt.hist(times)
plt.show()