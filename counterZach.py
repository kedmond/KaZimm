'''
Created May 15, 2009
Bernard Beckerman

This program retrieves a signal from the GPIB port using daq.py, finds the maximum
frequency of the signal using an FFT algorithm, and returns this value. The program
can also print and/or plot the data if the corresponding booleans are set to True.
---
Input: duration - time of data retrieval in seconds
Output: freq - maximum frequency of signal
'''

'''
Edited by Zachary Forbes, starting 20 Feb 2012
We wanted to be able to save the raw rotor frequencies, so I imported file naming and saving conventions directly from Zimm.py.
'''

import pylab
from daq import getData
from Zimm import initFile, getFilename, getHeader, writeToFile
from scipy.fftpack import fft, fftshift
import numpy
import matplotlib.pyplot as plt

def counter(duration, sampRate = 2048, printDataBool = False, plotDataBool = True, saveData = False, dataType = 'Rotor Frequencies'):
    data = retrieveData(duration, sampRate, saveData)     #Retrieve data from GPIB port, gives voltage as a function of sample number (sample number = time*sampRate) (see below)
    printData(printDataBool, data)              #Print data if printDataBool is set to True (see below)
    plotData(plotDataBool, data)                #Plot data if plotDataBool is set to True (see below)
    freq = analyzeData(data, duration)          #Computes frequency (see below)
    if saveData == True:
        initFile(dataType)                      #initializes file if we want to save the frequencies from this run
    return(freq)

def retrieveData(duration, sampRate, saveData):
    data = getData(int(sampRate), duration) #Takes data for duration seconds at sampRate samples per second (see daq.py)
    if saveData == True:
        writeToFile(data)                       #saves data if we want it
    return data

def printData(printDataBool, data):
    if printData:
        f = open('data.txt', 'w')
        for line in data:
            f.write(str(line)+'\n')
        f.close()
        
def analyzeData(data, duration):
    data = data - numpy.average(data)
    DATA = fftshift(fft(data))
    DATAsm = DATA*DATA.conj()
    m = max(DATAsm)
    N = len(DATAsm)
    for i in range (N/2, N):
        if DATAsm[i] == m: mi = i
    l = mi - N/2
    freq = (1.0*l/(duration))* (2.0 * numpy.pi / 100.0)*(180/numpy.pi)  #numerical factor 2pi radians per spoke, convert to degrees (180/pi)
    return freq

def plotData(plotDataBool, data):
    plt.ion()       # Turn on interactive mode.
    plt.hold(False) # Clear the plot before adding new data.
    plt.clf()
    plt.plot(data)
    plt.draw()
    plt.ioff()

'''
def retrieveDataOld(duration, sampRate):
    samples = duration*sampRate    
    n10s = int(pylab.floor(duration/10))
    remainingTime = int(duration%10)
    dataList = []
    if n10s >= 1:
        for i in range(n10s):
            dataList.append(getData(int(sampRate), 10))
    if remainingTime >= 1:
        dataList.append(getData(int(sampRate), remainingTime))
    data = []
    for i in range(n10s+1):
        data = pylab.append(data, dataList[i])
    return data

    count = 0
    su = 2.6 #upper hysteresis bound
    sl = 2.4 #lower hysteresis bound
    schmitt = 1
    i = 0
    countIndex = pylab.zeros(counts)
    sslc = pylab.zeros(counts) #samples since last count 

    
    return(data)
    while count < 1:
        i = i + 1
        if data[i] < sl:
            schmitt = 0
        if data[i] > su and schmitt == 0:
            #p = plt.axvline(x=i)
            count = count + 1
            schmitt = 1
        startSample = i
    count = 0
    for j in range (startSample, samples):
        sslc[count] = sslc[count] + 1 #increases variable sslc(samples since last count) by 1
        if data[j] < sl:
            schmitt = 0
        if data[j] > su and schmitt == 0:
            #p = plt.axvline(x=i)
            count = count + 1
            schmitt = 1
        endSample = j
        if count == counts:
            break
    goodCounts = 1
    for k in range (0, counts):
        if (1.5*pylab.average(sslc)) < sslc[k] or (.75*pylab.average(sslc)) > sslc[k]:
            print 'Count %d occurred at irregular interval.' % (k)
            goodCounts = 0
    samplesUsed = endSample - startSample
    timeElapsed = float(samplesUsed)/float(sampRate)
    Freq = float(counts)/timeElapsed
    if goodCounts == 0:
        #print 'bad counts'
        return 0
    else:
        return Freq
'''