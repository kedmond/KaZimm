'''
Created April 11th, 2014
Kazem Edmond
Based on KaZimm3

Modified by Kazem Edmond, 2011-

This program controls the Zimm-Crothers Viscometer and exports readings
of shear rate, shear stress, and viscosity to a text file.

Inputs:
startTemp - the temperature to start measurement
endTemp - the temperature to end measurement
incTemp - temperature difference between measurements
targMinSR - target minimum shear rate
targMaxSR - target maximum shear rate
numMotorFreqs - number of frequencies to be used in each temperature measurement
numMeas - number of measurements to take at each frequency
substance - string indicating substance name to use in Filename

'''

import bathControl
import counter as cnt
import peakDetect as pd

import serial
import glob
import numpy as np
import os
import sys

#import csv
import time
import re



def dataSet(temps, motorFreqs, durations=10.0, numMeas=10, pauseTime=1, sampRates=2048, 
            tempTime=120, motorTime=60, printDataBool=False, plotResults=False, sampleDescription='Water'):
    
    initialize()                   #Initialize ports, files, and constants (see below) 
    
    # save a description of the experiment to its directory
    descriptionHdr = 'duration (s), sampRate (cnts/s), motorFreq (1/s), pauseTime (s)'
    
    descriptionArray = np.column_stack((durations, sampRates, motorFreqs, np.repeat(pauseTime, durations.size)))
    
    directory = createDirectories()         # create a directory for this experiment's data
    
    # save a description of the experiment to its directory
    np.savetxt(directory + 'description-ViscosityTest.txt', descriptionArray, fmt='%d', delimiter=', ', header=sampleDescription+'\n'+descriptionHdr)    
    
    # these arrays need to be the same length
    if not motorFreqs.size == sampRates.size == durations.size:
        sys.exit('Error: motorFreqs, sampRates, and durations are not the same length.')
        
    for temp in temps:
        stirFreq = 2000.0                                    #Set stirring speed for bath
        setBath(temp, stirFreq, minStableTime=tempTime)     #Set temperature and wait for stability (see below)
        
        for i in range(motorFreqs.size):                    #Cycle through motorFreqs
            motorFreq = motorFreqs[i]
            print 'motorFreq is ' + str(motorFreq) + ' deg/sec.'
            
            sampRate = sampRates[i]
            
            print 'sampRate is ' + str(sampRate) + ' samples/sec.'
            duration = durations[i]
            # before each measurement, we reset the chains by stirring them up, and then letting them reform

            print 'Spinning rotor for ' + str(motorTime) + ' seconds.'             
            setMotorFreq(stirFreq, sleepTime=motorTime)     # Prevent sedimentation
            
            print 'Pausing rotor for ' + str(pauseTime) + ' seconds.'             
            setMotorFreq(0, sleepTime=pauseTime)            # Pause, perhaps to wait for chains to form 
            
            # Set motor and checks if motorFreq is too high
            breakBool = setMotorFreq(motorFreq, sleepTime=motorTime)
            # (high value may be changed in method setMotorFreq) (see below)
            if breakBool: break                             #Break if motorFreq is above threshold

            print 'Acquiring '+ str(numMeas) +' sets, each for a duration of ' + str(duration) + ' seconds.'
            
            # This will repeat the acquisition numMeas times,
            # This will repeat the acquisition numMeas times,
            # but then return all the data together in an array of stacked columns
            #Data retrieval (see below)
            rotorCounts = getData( numMeas, duration, sampRate ) 
            
            # calculate the std. error:
            # se = numpy.std(rotorFreqs) / sqrt(numpy.size(rotorFreqs))
            
            fname = 'rotorSpeeds-'+ str(motorFreq).zfill(4) +'.txt'
            
            # Data recording: temp. - motorFreq. - avg. rotor Freq. - std. err. - sampRate - duration
            writeToFile( rotorCounts, directory + fname, 
                        [duration, sampRate, motorFreq, motorTime, pauseTime] )
            
        print 'Measurements at ' + str(temp) + ' degrees C completed.'
        
    closeUp()


# Gather data while approaching a speed, and while slowing down again.
def rampSet(temps, motorFreqs, durations=10.0, numMeas=10, pauseTime=1, sampRates=2048, 
            tempTime=60, motorTime=60, printDataBool=False, plotResults=False, sampleDescription='Water'):
    
    initialize()                   #Initialize ports, files, and constants (see below) 
 
    # save a description of the experiment to its directory
    descriptionHdr = 'duration (s), sampRate (cnts/s), motorFreq (1/s), pauseTime (s)'
    
    descriptionArray = np.column_stack((durations, sampRates, motorFreqs, np.repeat(pauseTime, durations.size)))
    
    directory = createDirectories()         # create a directory for this experiment's data
    
    # save a description of the experiment to its directory
    np.savetxt(directory + 'description-RampTest.txt', descriptionArray, fmt='%d', delimiter=', ', header=sampleDescription+'\n'+descriptionHdr)    
    
    # these arrays need to be the same length
    if not motorFreqs.size == sampRates.size == durations.size:
        sys.exit('Error: motorFreqs, sampRates, and durations are not the same length.')
        
    for temp in temps:
        stirFreq = 2000.0                                   #Set stirring speed for bath
        setBath(temp, stirFreq, minStableTime=tempTime)     #Set temperature and wait for stability (see below)
        
        for i in range(motorFreqs.size):                    #Cycle through motorFreqs
            motorFreq = motorFreqs[i]
            sampRate = sampRates[i]
            duration = durations[i]
            
            print 'motorFreq is ' + str(motorFreq) + ' deg/sec, sampRate is ' + str(sampRate) + ' samples/sec, and duration is ' + str(duration) + ' seconds.'
            
            # before each measurement, we reset the chains by stirring them up, and then letting them reform

            rotorCountsUp = np.zeros([sampRate*duration, numMeas])
            rotorCountsDown = np.zeros([sampRate*duration, numMeas])
            
            for j in range(numMeas):
                
                # ***************** STIR UP SAMPLE *****************
                
                print 'Spinning rotor for ' + str(motorTime) + ' seconds.'             
                breakBool = setMotorFreq(stirFreq, sleepTime=motorTime)     # Prevent sedimentation
                if breakBool: break                             #Break if motorFreq is above threshold
                
                # ***************** REST SAMPLE *****************            
                
                print 'Pausing rotor for ' + str(pauseTime) + ' seconds.'             
                breakBool = setMotorFreq(0, sleepTime=pauseTime)            # Pause, perhaps to wait for chains to form 
                if breakBool: break                             #Break if motorFreq is above threshold
                
                # ***************** SHEAR SAMPLE FOR DATA *****************            
                
                # Must begin acquiring immediately after motor begins
                
                breakBool = setMotorFreq(motorFreq, sleepTime=0.0)
                
                rotorCountsUp[:, j] = cnt.retrieveData(duration, sampRate)
                
                rotorFreq = processData( rotorCountsUp[:, j], sampRate, duration )
                
                print 'Ramp up ' + str(j + 1) + ' of ' + str(numMeas) + ' taken successfully. [' + str(rotorFreq) + ' degs/sec]'
                
                # Wait for 15 seconds for steady shear
                time.sleep(15.0)
                
                # Stop the motor, and record how the sample slows
                breakBool = setMotorFreq(0.0, sleepTime=0.0)
                # (high value may be changed in method setMotorFreq) (see below)
                if breakBool: break                             #Break if motorFreq is above threshold
                
                rotorCountsDown[:, j] = cnt.retrieveData(duration, sampRate)
                
                rotorFreq = processData( rotorCountsDown[:, j], sampRate, duration )
                
                print 'Ramp down ' + str(j + 1) + ' of ' + str(numMeas) + ' taken successfully. [' + str(rotorFreq) + ' degs/sec]'

            fnameUp = 'rotorSpeedsUp-'+ str(motorFreq).zfill(4) +'.txt'
            fnameDown = 'rotorSpeedsDown-'+ str(motorFreq).zfill(4) +'.txt'
            
            # Data recording: temp. - motorFreq. - avg. rotor Freq. - std. err. - sampRate - duration
            writeToFile( rotorCountsUp, directory + fnameUp, 
                        [duration, sampRate, motorFreq, motorTime, pauseTime] )
                        
            # Data recording: temp. - motorFreq. - avg. rotor Freq. - std. err. - sampRate - duration
            writeToFile( rotorCountsDown, directory + fnameDown, 
                        [duration, sampRate, motorFreq, motorTime, pauseTime] )
            
        print 'Measurements at ' + str(temp) + ' degrees C completed.'
        
    closeUp()


# Gather data without stopping the motor in between sets
def hystSet(temps, motorFreqs, durations=10.0, sampRates=2048,
            tempTime=60, motorTime=15, sampleDescription='Water'):
    
    revFreqs = motorFreqs[(motorFreqs.size-2)::-1]
    revRates = sampRates[(motorFreqs.size-2)::-1]
    revDurs = durations[(motorFreqs.size-2)::-1]
    
    initialize()                   #Initialize ports, files, and constants (see below) 
 
    # save a description of the experiment to its directory
    descriptionHdr = 'duration (s), sampRate (cnts/s), motorFreq (1/s), pauseTime (s)'
    
    descriptionArray = np.column_stack((durations, sampRates, motorFreqs, np.repeat(0, durations.size)))
    
    directory = createDirectories()         # create a directory for this experiment's data
    
    # save a description of the experiment to its directory
    np.savetxt(directory + 'description-HysterisisTest.txt', descriptionArray, fmt='%d', delimiter=', ', header=sampleDescription+'\n'+descriptionHdr)    
    
    # these arrays need to be the same length
    if not motorFreqs.size == sampRates.size == durations.size:
        sys.exit('Error: motorFreqs, sampRates, and durations are not the same length.')
        
    for temp in temps:
        stirFreq = 2000.0                                   #Set stirring speed for bath
        setBath(temp, stirFreq, minStableTime=tempTime)     #Set temperature and wait for stability (see below)
        
        for i in range(motorFreqs.size):                    #Cycle through motorFreqs
            motorFreq = motorFreqs[i]
            sampRate = sampRates[i]
            duration = durations[i]
            
            print 'motorFreq is ' + str(motorFreq) + ' deg/sec, sampRate is ' + str(sampRate) + ' samples/sec, and duration is ' + str(duration) + ' seconds.'
                        
            # ***************** SHEAR SAMPLE FOR DATA *****************            
            
            # We want to measure steady shear, so wait for motor and sample to achieve speed
            breakBool = setMotorFreq(motorFreq, sleepTime=motorTime)
            if breakBool: break                             #Break if motorFreq is above threshold
            
            rotorCounts = cnt.retrieveData(duration, sampRate)
            
            rotorFreq = processData( rotorCounts, sampRate, duration )
            
            print 'Data for '+ str(motorFreq) +' deg/sec taken successfully. [rotorFreq: ' + str(rotorFreq) + ' degs/sec]'
            
            # Wait for 15 seconds for steady shear
            #time.sleep(15.0)
            
            fname = 'rotorSpeedsUp-'+ str(motorFreq).zfill(4) +'.txt'
            
            # Data recording: temp. - motorFreq. - avg. rotor Freq. - std. err. - sampRate - duration
            writeToFile( rotorCounts, directory + fname, 
                        [duration, sampRate, motorFreq, motorTime, 0] )
        
        
        for i in range(revFreqs.size):                    #Cycle through motorFreqs
            motorFreq = revFreqs[i]
            sampRate = revRates[i]
            duration = revDurs[i]
            
            print 'motorFreq is ' + str(motorFreq) + ' deg/sec, sampRate is ' + str(sampRate) + ' samples/sec, and duration is ' + str(duration) + ' seconds.'
                        
            # ***************** SHEAR SAMPLE FOR DATA *****************            
            
            # We want to measure steady shear, so wait for motor and sample to achieve speed
            breakBool = setMotorFreq(motorFreq, sleepTime=motorTime)
            if breakBool: break                             #Break if motorFreq is above threshold
            
            rotorCounts = cnt.retrieveData(duration, sampRate)
            
            rotorFreq = processData( rotorCounts, sampRate, duration )
            
            print 'Data for '+ str(motorFreq) +' deg/sec taken successfully. [rotorFreq: ' + str(rotorFreq) + ' degs/sec]'
            
            # Wait for 15 seconds for steady shear
            #time.sleep(15.0)
            
            fname = 'rotorSpeedsDown-'+ str(motorFreq).zfill(4) +'.txt'
            fname = 'rotorSpeedsDown-'+ str(motorFreq).zfill(4) +'.txt'
            
            # Data recording: temp. - motorFreq. - avg. rotor Freq. - std. err. - sampRate - duration
            writeToFile( rotorCounts, directory + fname, 
                        [duration, sampRate, motorFreq, motorTime, 0] )
            
        print 'Measurements at ' + str(temp) + ' degrees C completed.'
        
    closeUp()


def writeToFile(rotorCounts, fileName, parameters):
    # a list of pertinent parameters used in the experiment
    duration = parameters[0]
    sampRate = parameters[1]
    motorFreq = parameters[2]
    motorTime = parameters[3]
    pauseTime = parameters[4]
    
    rotorFreqs = processData( rotorCounts, sampRate, duration ) #Data processing (see below)
    
    # put the rotor speeds as the first row of each column of data
    data = np.insert(rotorCounts, 0, rotorFreqs, axis=0)
    
    np.savetxt( fileName, data, delimiter=', ', header='duration: ' + 
    str(duration) + ', sampRate: ' + str(sampRate) + ', motorFreq: ' + 
    str(motorFreq) + ', motorTime: ' + str(motorTime) + ', pauseTime: ' + 
    str(pauseTime) )


# Create directory structure for data from experiments
def createDirectories(rootDir=os.getcwd()):
    
    # make sure it ends with a slash
    if not (rootDir[-1] == os.sep):
        rootDir += os.sep
    
    dateDir = time.strftime('%d%b%Y') + os.sep
    expStr = 'Experiment'
    expNum = 1
    
    # Look for other "Experiment" directories, and create the next one 
    
    if not os.path.exists(rootDir + dateDir + expStr + str(expNum)):
        os.makedirs(rootDir + dateDir + expStr + str(expNum))
        # going forward, this is where to store data
        directory = rootDir + dateDir + expStr + str(expNum)
    else:
        # get a list of directories that start with 'Experiment', count them, and enumerate
        fns = glob.glob(rootDir + dateDir + expStr + '*')
        fns.sort()
        
        # figure out what number this new directory should be
        expNum = int(re.findall(r'\d+', fns[-1][-11:])[0]) + 1
        
        # Make a directory in the correct spot, with a name
        os.makedirs(rootDir + dateDir + expStr + str(expNum))
        
        # going forward, this is where to store data
        directory = rootDir + dateDir + expStr + str(expNum)
    
    # End the directory string with the appropriate separator
    directory += os.sep
    
    return directory



def loadData(directory):
    # assume each file has the same number of columns
    # calculate the frequency of each column and its std. dev., for error bars
    # record the motorFreqs from the description file
        
    if directory[-1] != os.sep: directory += os.sep
    
    fileNames = glob.glob(directory+'rotorSpeeds*')
    descriptionFile = glob.glob(directory + 'description*')[0]

    # extract the motorFreqs from the first column    
    motorFreqs = np.genfromtxt(descriptionFile, delimiter=',')[:, 2]
    durations = np.genfromtxt(descriptionFile, delimiter=',')[:, 0]
    
    # two columns: motor speeds and rotor speeds
    data = np.zeros((motorFreqs.size, 2))
    data[:, 0] = motorFreqs
    
    for file in fileNames:
        data = np.genfromtxt(file, skiprows=2, delimiter=',')
        
        ncols = data.shape[1]
        for i in ncols:
            means = np.mean(data[:, i])
            
            
            
            


def processData(rotorCounts, sampRate, duration):
    
    # The number of columns tells us the number of data sets
    if rotorCounts.ndim == 1: 
        rotorFreqs = cnt.analyzeData(rotorCounts, duration)
    else:        
        numMeas = rotorCounts[0, :].size
        
        # store the rotor frequences here
        rotorFreqs = np.zeros(numMeas)
        
        for i in xrange(numMeas):
            rotorFreqs[i] = cnt.analyzeData(rotorCounts[:, i], duration)
        
    return rotorFreqs


# Computes data (visc, ss, sr) from angular frequencies of inner tube 
# and angular frequency of motor, concatenates data into string, and
# returns data (called by main)
def processDataOld(temp, motorFreq, rotorFreqs):
    # Convert motor speed from degrees per second to radians per second
    motorFreqRads = motorFreq * np.pi / 180.0    #only important if using rad/s for rotorFreqs!
    
    rotorFreq = np.mean(rotorFreqs)                         # Average all measurements
    
    # This array of strings contains the data.
    se = np.std(rotorFreqs) / np.sqrt(np.size(rotorFreqs))
    data = str(temp)+','+str(motorFreqRads)+','+str(rotorFreq)+','+str(se)
    
    return data                                            # Returns data string


def initialize():    #Method initializes motor, bath, and constants to be used in calculation
    initMotor()                 #Initialize motor
    initBath()                  #Initialize thermal bath
    setConsts()                 #Set constants to be used in calculation (must be edited if tubes are changed)
        

def initMotor():#Initializes motor for use in experiment (called by initialize)
    global Motor                #Initialize global variable Motor for use as name of the motor port
    Motor = serial.Serial('COM6')#Initialize serial port and assign to it the global variable Motor
#    Motor.open()                #Open communication with the port labeled Motor
    Motor.write('EN\r\n')       #Enable motor using command explained in Manual.doc
        
    
def initBath():#Initializes bath for use in experiment (called by initialize)
    global TBath                #Initialize global variable TBath for use as name of the Temperature Bath
    TBath = bathControl.RW()    #Initialize object TBath (see bathControl.py)
    

def setConsts():#Sets constants for use in calculation (called by initialize)
    global GeoConst             #Initialize global GeoConst - geometrical constant for use in calculation
    global CalConst             #Initialize global CalConst - calibration constant for use in calculation
    ri = 12.0                   #Outer radius of the inner tube in mm
    ro = 12.95                  #Inner radius of the outer tube in mm
    
    GeoConst = 8.0 * np.pi * ro**2.0 * ri**2.0 * np.log10(ro / ri) / (ro**2.0 - ri**2.0)**2.0 #Compute GeoConst
    CalConst = 2.0*5.965252117              #Set CalConst - must be recalibrated if tubes are switched


#Sets bath temperature and waits for equilibration (called by main)
def setBath(temp, stirFreq, minStableTime=120.0):
    
    TBath.setSet(temp)          # Set setpoint temperature (see bath.py)
    setMotorFreq(stirFreq, sleepTime=2.0)  # Set a low frequency shear while temperature equilibrates
    print 'Waiting for temp to stabilize at ' + str(temp) + ' degrees C.'
    unstable = True
    #minStableTime = 120.0       #Minimum time required for the bath to be at steady stayed before proceeding
    while unstable:
        t = time.time()
        unstable = False
        while abs(time.time() - t) < minStableTime: # Waits for stability for minStableTime
            if abs(TBath.extTemp() - temp) > .01:# Checks for temperature equal to the setpoint
                unstable = True
                break           #Restarts count if the bath temperature is not equal to setpoint at any point during the count
    print 'Temp = ' + str(temp) + ' degrees C. Stable setpoint reached.'
    # If minStableTime is reached, a message is printed and the program moves on



# Goes through all shear rates and samples each one numMeas times
# Returns a stack of columns 
def getData(numMeas, duration, sampRate):
    
    rotorCounts = np.zeros([sampRate*duration, numMeas])
    
    for j in range(0, numMeas):
        rotorCounts[:, j] = cnt.retrieveData(duration, sampRate)
        
        rotorFreq = processData( rotorCounts[:, j], sampRate, duration )
        
        print 'count ' + str(j + 1) + ' of ' + str(numMeas) + ' taken successfully. [' + str(rotorFreq) + ' degs/sec]'
        # Kazem, 3.16.2012 - I added the reporting of data for convenience.
    
    return rotorCounts


# Goes through all shear rates and samples each one numMeas times; 
# returns angular frequencies of inner tube (called by main)
def getDataOld(numMeas, countTime=10.0, sampRate=2048, printDataBool=False):
    
    rotorFreqs = np.zeros(numMeas)          # Create array in which to store rotorFreqs
    
    for j in range (0, numMeas):            # Take numMeas readings
        while rotorFreqs[j] == 0:           # It is possible to have an unsuccessful count.
                                            # This checks if any previous counts have been 
                                            # successful and allows another count if necessary

            #edit 24/1/2013 by Zach
            #save counter data into new file            
            count = cnt.counter(countTime, plotDataBool=False, sampRate=sampRate, printDataBool=printDataBool)            
            rotorFreqs[j] = count #* math.pi / 180
            #Angular frequency measurement of inner tube (numerical factor 2pi/100 radians per spoke)
            #writeToFile(count)#, f = 'counter-'+Filename)
            
        print 'count ' + str(j + 1) + ' of ' + str(numMeas)
        + ' taken successfully. [' + str(rotorFreqs[j]) + ' degs/sec]'
        # Kazem, 3.16.2012 - I added the reporting of data for convenience.
    return rotorFreqs         #Return angular frequencies of inner tube
            

def closeUp():#Ends Experiment (called by main)
    print 'Experiment Complete.'
    #endTime = time.strftime("%y-%m-%d-%H-%M-%S", time.localtime())#System time
    #writeToFile(endTime)                                    #Prints end time to file
    TBath.close()                                           #Closes temperature bath
    Motor.write('DI\r\n')                                   #Disables motor (see Manual.doc for details)
    Motor.close()                                           #Closes motor
    
# Checks to see if motorFreq is small enough to be safe, sets the motor frequency,
# and waits for rotor to speed up (called by main, setBath)
def setMotorFreq(motorFreq, sleepTime=60):  
    skipBool = False
    
    if motorFreq > 2000.0:                  # Check motorFreq is safe
        print str(motorFreq) + ' deg/s too large. Measurement(s) skipped.'
        skipBool = True                     # Set skipBool to True: tell program to skip data analysis for this measurement
    else:
        Motor.write('FR ' + str(motorFreq) + '\r\n')
        # If motorFreq is safe, set motor to freerun at motorFreq (see Manual.doc for details on setting motor speed)
        print 'Motor frequency = ' + str(motorFreq) + ' degrees per second (' + str(sleepTime) + ' second pause).'
        time.sleep(sleepTime)               # Wait for rotor to speed up
    return skipBool                         # Returns skipBool that tells main to skip data analysis
    
