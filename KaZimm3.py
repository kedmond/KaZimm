'''
Created May 28th, 2009
Bernard Beckerman

Modified by Kazem Edmond, 2011-

25.Feb.2014 - Kaz added parameter to pause rotor to allow colloidal chains to reform before shearing
28.Feb.2014 - Kaz lets you define a different count time for each set of data

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

To do: let the user enter a list of count rates and durations.

'''

from __future__ import division
#from numarray import *
from pylab import *
#import pylab
import time
import math
import counter
import serial
import bathControl
import numpy
import kplotData3 as kplt
import os

# read/write to Excel spreadsheet 
#import xlwt
#import xlrd


def main(startTemp=40, endTemp=17, incTemp=-0.5, targMinSR = 1.0, targMaxSR = 21.0, numMotorFreqs = 11, numMeas=5, substance = 'Water'):
    initialize(substance)                                   #Initialize ports, files, and constants (see below) 
    temps = range(startTemp, endTemp + incTemp, incTemp)   #Set temperatures to be used
    for temp in temps:                                      #Cycle through temps
        stirFreq = 200.0                                    #Set stirring speed for bath
        setBath(temp, stirFreq)                             #Set temperature and wait for stability (see below)
        motorFreqs = motorCalibration(targMinSR, targMaxSR, numMotorFreqs, stirFreq) #Set motorFreqs to be used (see below)
        for motorFreq in motorFreqs:                        #Cycle through motorFreqs
            breakBool = setMotorFreq(motorFreq)             #Set motor and checks if motorFreq is too high (high value may be changed in method setMotorFreq) (see below)
            if breakBool: break                             #Break if motorFreq is above threshold
            rotorFreqs = getData(numMeas)                   #Data retrieval (see below)
            data = processData(temp, motorFreq, rotorFreqs) #Data processing (see below)
            writeToFile(data)                               #Data recording (see below)
        print 'Measurements at ' + str(temp) + ' degrees C completed. endTemp = ' + str(endTemp) + '.'
    closeUp()                                               #Closes up (see below)


# 2.15.2012, Kazem: Added this to just do a set of rates at one temperature.
def singleSet(temps, motorFreqs, autoCount=False, countTime=10.0, numMeas=10, pauseTime=1, sampRate=2048, tempTime=120, motorTime=60, printDataBool=False, plotResults = True, substance='Water'):
    initialize(substance)                               #Initialize ports, files, and constants (see below) 
    
    for temp in temps:
        stirFreq = 400.0                                    #Set stirring speed for bath
        setBath(temp, stirFreq, minStableTime=tempTime)     #Set temperature and wait for stability (see below)
        
        for motorFreq in motorFreqs:                        #Cycle through motorFreqs
            print 'motorFreq is ' + str(motorFreq) + ' deg/sec.'
            # 28.Feb.2014, Kaz: minimize sampRate to avoid data overflow for longer runs at slow shear:
            rates = 2**np.array(range(4, 12))
            
            minSampRate = (20.0 * motorFreq/3.6)
            # rational is that you need at least 15 counts per tick, which I determined experimentally
            # with 100 ticks per rotation, and 3.6 degrees per tick
            # 15 counts/tick * (motorRate deg./sec. / 3.6 deg./tick) = minCountRate counts/second
            
            sampRate = rates[np.abs(rates-minSampRate).argmin()]
            if (sampRate-minSampRate) < 0 and sampRate != rates[-1]:
                sampRate = rates[np.abs(rates-minSampRate).argmin()+1]
         	
            print 'sampRate is ' + str(sampRate) + ' samples/sec.'
            
            # before each measurement, we reset the chains by stirring them up, and then letting them reform
            print 'Spinning rotor for ' + str(motorTime) + ' seconds.'             
            setMotorFreq(stirFreq, sleepTime=motorTime)    # Prevent sedimentation
            
            print 'Pausing rotor for ' + str(pauseTime) + ' seconds.'             
            setMotorFreq(0, sleepTime=pauseTime)    # Pause, perhaps to wait for chains to form 
            
            breakBool = setMotorFreq(motorFreq, sleepTime=motorTime) #Set motor and checks if motorFreq is too high (high value may be changed in method setMotorFreq) (see below)
            if breakBool: break                             #Break if motorFreq is above threshold
            
            # autoCount sets the count time to however long it takes for the motor to complete to full rotations
            if autoCount == True:
                cTime = np.ceil(720./motorFreq)
                print 'Counting for ' + str(cTime) + ' seconds.'
                rotorFreqs = getData(numMeas, countTime=cTime, sampRate=sampRate, printDataBool=printDataBool)     #Data retrieval (see below)
            else:
                rotorFreqs = getData(numMeas, countTime=countTime, sampRate=sampRate, printDataBool=printDataBool) #Data retrieval (see below)
                
            plt.savefig('Figs/'+time.strftime("%y-%m-%d-%H-%M-%S", time.localtime())+'-'+str(motorFreq)+'.png')
            data = processData(temp, motorFreq, rotorFreqs) #Data processing (see below)
            
            # calculate the std. error:
            #se = numpy.std(rotorFreqs) / sqrt(numpy.size(rotorFreqs))
            
            # Data recording: temp. - motorFreq. - avg. rotor Freq. - std. err. - sampRate - duration
            writeToFile( data )
            
        print 'Measurements at ' + str(temp) + ' degrees C completed.'
        
    closeUp()   
    
    if plotResults and len(data) > 1:
        kplt.plotData(os.getcwd(), Filename, writeFile=True, fit=False)

# 28.March.2014, Kaz: I want to define different sample rates for each frequency.
def dataSet(temps, motorFreqs, durations=10.0, numMeas=10, pauseTime=1, sampRates=2048, tempTime=120, motorTime=60, printDataBool=False, plotResults = True, substance='Water' ):
    initialize(substance)                               #Initialize ports, files, and constants (see below) 
    
    for temp in temps:
        stirFreq = 400.0                                    #Set stirring speed for bath
        setBath(temp, stirFreq, minStableTime=tempTime)     #Set temperature and wait for stability (see below)
        
        for i in range(motorFreqs.size):                        #Cycle through motorFreqs
            motorFreq = motorFreqs[i]
            print 'motorFreq is ' + str(motorFreq) + ' deg/sec.'
            
            sampRate = sampRates[i]
            
            print 'sampRate is ' + str(sampRate) + ' samples/sec.'
            duration = durations[i]
            # before each measurement, we reset the chains by stirring them up, and then letting them reform

            print 'Spinning rotor for ' + str(motorTime) + ' seconds.'             
            setMotorFreq(stirFreq, sleepTime=motorTime)    # Prevent sedimentation
            
            print 'Pausing rotor for ' + str(pauseTime) + ' seconds.'             
            setMotorFreq(0, sleepTime=pauseTime)    # Pause, perhaps to wait for chains to form 
            
            breakBool = setMotorFreq(motorFreq, sleepTime=motorTime) #Set motor and checks if motorFreq is too high (high value may be changed in method setMotorFreq) (see below)
            if breakBool: break                             #Break if motorFreq is above threshold
            
            rotorFreqs = getData(numMeas, countTime=duration, sampRate=sampRate, printDataBool=printDataBool) #Data retrieval (see below)

            plt.savefig('Figs/'+time.strftime("%y-%m-%d-%H-%M-%S", time.localtime())+'-'+str(motorFreq)+'.png')
            
            data = processData(temp, motorFreq, rotorFreqs) #Data processing (see below)
            
            # calculate the std. error:
            #se = numpy.std(rotorFreqs) / sqrt(numpy.size(rotorFreqs))
            
            # Data recording: temp. - motorFreq. - avg. rotor Freq. - std. err. - sampRate - duration
            writeToFile( data )
            
        print 'Measurements at ' + str(temp) + ' degrees C completed.'
        
    closeUp()   
    
    if plotResults and len(data) > 1:
        kplt.plotData(os.getcwd(), Filename, writeFile=True, fit=False)


def initialize(substance):    #Method initializes motor, bath, file, and constants to be used in calculation
    initMotor()                 #Initialize motor
    initBath()                  #Initialize thermal bath
    initFile(substance)         #Initialize and name file
    setConsts()                 #Set constants to be used in calculation (must be edited if tubes are changed)


def initMotor():#Initializes motor for use in experiment (called by initialize)
    global Motor                #Initialize global variable Motor for use as name of the motor port
    Motor = serial.Serial('COM6')#Initialize serial port and assign to it the global variable Motor
#    Motor.open()                #Open communication with the port labeled Motor
    Motor.write('EN\r\n')       #Enable motor using command explained in Manual.doc
    
def initBath():#Initializes bath for use in experiment (called by initialize)
    global TBath                #Initialize global variable TBath for use as name of the Temperature Bath
    TBath = bathControl.RW()    #Initialize object TBath (see bathControl.py)


def initFile(substance):#Initializes file and writes title and column names for use in data recording (called by initialize)
    global Filename             #Initialize global variable Filename
    Filename = getFilename(substance)#Create Filename
    header = getHeader(substance)#Create header  
    writeToFile(header)         #Write header to file


def getFilename(substance):#Creates Filename with substance and time (called by initFile)
    startTime = time.strftime("%y-%m-%d-%H-%M-%S", time.localtime()) #create YY-MM-DD-HH-MM-SS timestring
    filename = '%s_%s.txt' % (startTime, substance) #Set filename
    return filename

def getHeader(substance):#Creates header (title and column names) using substance (called by initFile)
    line1 = 'Shear Rate, Shear Stress, and Viscosity Measurements with the Zimm-Crothers Viscometer: %s\n' % (substance)
    line2 = 'Temp,MotorFrequency,RotorFrequency'
    header = line1 + line2
    return header

#def writeToFile(line, f = fname):# Writes a string to file (called by main, initFile, closeUp)
    #edit 24/1/2013 by Zach
    #added f parameter to allow for writing a new file with frequencies gathered in getData   
def writeToFile(line):# Writes a string to file (called by main, initFile, closeUp)
    file = open(Filename, 'a')      #Open file to append
#    file = open(f, 'a')            #Open file to append
    file.write(line + '\n')         #Write line + escape char
    file.close()                    #Close file so it can be accessed during data runs

def setConsts():#Sets constants for use in calculation (called by initialize)
    global GeoConst             #Initialize global GeoConst - geometrical constant for use in calculation
    global CalConst             #Initialize global CalConst - calibration constant for use in calculation
    ri = 12.0                   #Outer radius of the inner tube in mm
    ro = 12.95                  #Inner radius of the outer tube in mm
    GeoConst = 8.0 * math.pi * ro**2.0 * ri**2.0 * log(ro / ri) / (ro**2.0 - ri**2.0)**2.0 #Compute GeoConst
    CalConst = 2.0*5.965252117              #Set CalConst - must be recalibrated if tubes are switched

def setBath(temp, stirFreq, minStableTime=120.0):#Sets bath temperature and waits for equilibration (called by main)
    TBath.setSet(temp)          # Set setpoint temperature (see bath.py)
    setMotorFreq(stirFreq, sleepTime=2.0)  # Set a low frequency shear while temperature equilibrates
    print 'Waiting for temp to stabilize at ' + str(temp) + ' degrees C.'
    unstable = True
    #minStableTime = 120.0       #Minimum time required for the bath to be at steady stayed before proceeding
    while unstable:
        t = time.time()
        unstable = False
        while abs(time.time() - t) < minStableTime: #Waits for stability for minStableTime
            if abs(TBath.extTemp() - temp) > .01:#Checks for temperature equal to the setpoint
                unstable = True
                break           #Restarts count if the bath temperature is not equal to setpoint at any point during the count
    print 'Temp = ' + str(temp) + ' degrees C. Stable setpoint reached.'#If minStableTime is reached, a message is printed and the program moves on

# Calibrates the viscometer so that the minimum shear rate is close to 
# targMin, the maximum shear rate is close to maxSR, and the others are 
# equally spaced in between (called by main)
def motorCalibration(targMin, targMax, numMotorFreqs, motorCalFreq):
    rotorCalFreq = counter.counter(8, plotDataBool=False)* math.pi / 180#Calibration Measurement (numerical factor 2pi/100 radians per spoke)
    calShearRate = GeoConst * rotorCalFreq                  #Calculate shear rate from calibration run
    motorFreqMin = motorCalFreq * targMin / calShearRate    #Calculate Minimum Motor Freq (assume linear relation between motor speed and shear rate)
    motorFreqMax = motorCalFreq * targMax / calShearRate    #Calculate Maximum Motor Freq (assume linear relation between motor speed and shear rate)
    df = (motorFreqMax-motorFreqMin)/(numMotorFreqs-1)      #Use number of freqs to be used and min and max to determine difference between consecutive motor frequencies
    motorFreqs = arange(motorFreqMin, motorFreqMax+df, df)  #Create 1-d array of motor freqs to be used
    return motorFreqs                                       #Return motorFreqs

def calculateViscosity(temp, motorFreq, rotorFreqs):
    # Calculate viscosity from angular frequencies of inner tube and angular frequency of motor.
    # Returns value for viscosity.
    motorFreqRadians = motorFreq * math.pi / 180.0          #Convert motor speed from degrees per second to radians per second
    rotorFreq = average(rotorFreqs)                         #Average all measurements
    v = CalConst *(motorFreqRadians - rotorFreq) / (rotorFreq)#Calculate viscosity
    return v

# Computes data (visc, ss, sr) from angular frequencies of inner tube 
# and angular frequency of motor, concatenates data into string, and
# returns data (called by main)
def processData(temp, motorFreq, rotorFreqs):
    motorFreqRads = motorFreq * math.pi / 180.0          # Convert motor speed from degrees per second to radians per second
    #only important if using rad/s for rotorFreqs!
    rotorFreq = average(rotorFreqs)                         # Average all measurements
    
    # This array of strings contains the data.
    se = numpy.std(rotorFreqs) / sqrt(numpy.size(rotorFreqs))
    data = str(temp)+','+str(motorFreqRads)+','+str(rotorFreq)+','+str(se)
    
    return data                                            # Returns data string

# Goes through all shear rates and samples each one numMeas times; 
# returns angular frequencies of inner tube (called by main)
def getData(numMeas, countTime=10.0, sampRate=2048, printDataBool=False):
    rotorFreqs = zeros(numMeas)		               # Create array in which to store rotorFreqs
    
    for j in range (0, numMeas):                     # Take numMeas readings
        while rotorFreqs[j] == 0:                    # It is possible to have an unsuccessful count. This checks if any previous counts have been successful and allows another count if necessary
            #edit 24/1/2013 by Zach
            #save counter data into new file            
            count = counter.counter(countTime, plotDataBool=False, sampRate=sampRate, printDataBool=printDataBool)            
            rotorFreqs[j] = count#* math.pi / 180    #Angular frequency measurement of inner tube (numerical factor 2pi/100 radians per spoke)
            #writeToFile(count)#, f = 'counter-'+Filename)

        print 'count ' + str(j + 1) + ' of ' + str(numMeas) + ' taken successfully. [' + str(rotorFreqs[j]) + ' degs/sec]'
        # Kazem, 3.16.2012 - I added the reporting of data for convenience.
    return rotorFreqs         #Return angular frequencies of inner tube
            

def closeUp():#Ends Experiment (called by main)
    print 'Experiment Complete.'
    #endTime = time.strftime("%y-%m-%d-%H-%M-%S", time.localtime())#System time
    #writeToFile(endTime)                                    #Prints end time to file
    TBath.close()                                           #Closes temperature bath
    Motor.write('DI\r\n')                                   #Disables motor (see Manual.doc for details)
    Motor.close()                                           #Closes motor
    

def setMotorFreq(motorFreq, sleepTime=60):#Checks to see if motorFreq is small enough to be safe, sets the motor frequency, and waits for rotor to speed up (called by main, setBath)
    skipBool = False
    if motorFreq > 2000.0:                                  # Check motorFreq is safe
        print str(motorFreq) + ' deg/s too large. Measurement(s) skipped.'
        skipBool = True                                     # Set skipBool to True: tell program to skip data analysis for this measurement
    else:
        Motor.write('FR ' + str(motorFreq) + '\r\n')        # If motorFreq is safe, set motor to freerun at motorFreq (see Manual.doc for details on setting motor speed)
        print 'Motor frequency = ' + str(motorFreq) + ' degrees per second.'
        time.sleep(sleepTime)                               # Wait for rotor to speed up
    return skipBool                                         # Returns skipBool that tells main to skip data analysis
