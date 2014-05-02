'''
Created May 15, 2009
Bernard Beckerman

This program retrieves a signal from the GPIB port using daq.py, finds the maximum
frequency of the signal using an FFT algorithm, and returns this value. The program
can also print and/or plot the data if the corresponding booleans are set to True.
---
Input: duration - time of data retrieval in seconds
Output: freq - maximum frequency of signal

History:
    6/Nov/2013 - Kazem: I edited the daq.py code so that the software waits for 
    the full duration to acquire data.  Before, it was arbitrarily set to 10 seconds.
    I am now modifying the code to increase duration for low shear rates.
    
    20/Mar/2014 - Kaz: I changed the printData routine to include details about
    the aquisition, and it now uses numpy.savetxt()
'''
#import pylab
from time import gmtime, strftime
from daq import getData
from scipy.fftpack import fft, fftshift
import numpy as np
import matplotlib.pyplot as plt
import csv

def counter(duration, sampRate = 2048, printDataBool = False, plotDataBool = True):
    # Retrieve data from GPIB port, gives voltage as a function of sample number
    # (sample number = time*sampRate) (see below)
    data = retrieveData(duration, sampRate)
    
    if printDataBool:
        dataOut = np.column_stack( (np.array(range(0, data.size), dtype=int), data))
        
        printData(dataOut, duration, sampRate)  #Print data if printDataBool is set to True (see below)
        
    if plotDataBool:     
        plotData(data, sampRate)                #Plot data if plotDataBool is set to True (see below)
        
    freq = analyzeData(data, duration)          #Computes frequency (see below)
    
    return(freq)


def retrieveData(duration, sampRate):
# Kazem - 6/Nov/2013 - If shear rate is really low, duration needs to be increased
# Kaz - 28.Fev.2014 - long durations require lower count rates.  will fill up memory.

    data = getData(int(sampRate), duration) #Takes data for duration seconds at sampRate samples per second (see daq.py)
    return data


def printData(data, duration, sampRate):
    f = open('data-'+strftime("%Y-%m-%d_%H.%M", gmtime())+'.txt', 'w')
    
    f.write('duration: '+str(duration)+', sample rate: '+str(sampRate)+'\n')
  
    np.savetxt(f, data, delimiter=',', newline='\n', fmt='%i, %f')  

    f.close()
    
    
def analyzeData(data, duration):
    data = data - np.average(data)
    # Kaz, 18.Mar.2014
    # The fftshift is unnecessary.
    # I should use the fftpack.fftfreq routine instead
    
    DATA = fftshift(fft(data))
    DATAsm = DATA*DATA.conj()
    m = max(DATAsm)
    N = len(DATAsm)
    for i in range(N/2, N):
        if DATAsm[i] == m: mi = i
    l = mi - N/2
    freq = (1.0*l/(duration))* (2.0 * np.pi / 100.0)*(180/np.pi)
    
    #numerical factor 2pi radians per spoke, convert to degrees (180/pi)
    
    return freq
    
    
def analyzeFile(filename):
    
    # load up the data
    data = np.genfromtxt(filename, delimiter=',', skip_header=1)
    
    # the first row is just the row numbers
    data = data[:, 1]    
    
    # I don't think this is necessary...
    data = data - np.average(data)
    
    # get the duration from the file's header
    r = csv.reader(open(filename))
    line1 = r.next()
    values = [int(s) for s in (line1[0]+line1[1]).split() if s.isdigit()]
    duration = values[0]
    
    # Kaz, 18.Mar.2014
    # The fftshift is unnecessary.
    # I should use the fftpack.fftfreq routine instead
    
    DATA = fftshift(fft(data))
    DATAsm = DATA*DATA.conj()
    m = max(DATAsm)
    N = len(DATAsm)
    for i in range(N/2, N):
        if DATAsm[i] == m: mi = i
    l = mi - N/2
    freq = (1.0*l/(duration))* (2.0 * np.pi / 100.0)*(180/np.pi)
    
    #numerical factor 2pi radians per spoke, convert to degrees (180/pi)
    
    return freq

def plotData(data, sampRate):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(data)
    
    ax.set_xlim([0, data.size])
    ax.set_ylim([0, 1.1*data.max()])
    # 28.Feb.2014, Kaz: change x-axis to seconds
    labels = ax.get_xticks()
    labels /= sampRate
    ax.set_xticklabels(labels)

    ax.set_xlabel('duration (s)')
    
    fig.subplots_adjust(left=0.025, right=0.99, top=0.99, bottom=0.15)
    
    ax.grid(True)    
    
    fig.canvas.draw()
    
#def plotData(data, sampRate):
#    plt.ion()       # Turn on interactive mode.
#    #plt.hold(False) # Clear the plot before adding new data.
#    plt.clf()
#    plt.plot(data)
#    # 2.15.2012, Kazem: Added ylim/xlim for nicer looking plots.
#    #plt.ylim([0, 1.1*data.max()]) 
#    plt.ylim([0, 11.0])
#    plt.xlim([0, data.size])
#    # 28.Feb.2014, Kaz: change x-axis to seconds
#    plt.xticks(range(0, int(data.size+1), int(sampRate)), range(int(data.size/sampRate+1)))
#    plt.xlabel('duration (seconds)')
#    
#    plt.subplots_adjust(left=0.025, right=0.99, top=0.99, bottom=0.15)
#    
#    plt.grid(True)    
#    
#    plt.draw()
#    plt.ioff()

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