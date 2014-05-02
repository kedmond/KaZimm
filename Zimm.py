'''
Created May 28th, 2009
Bernard Beckerman

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

#from numarray import *
from pylab import *
import time
import math
import counter
import serial
import bathControl

def main(startTemp=40, endTemp=17, incTemp=-0.5, targMinSR = 1.0, targMaxSR = 21.0, numMotorFreqs = 11, numMeas=5, substance = 'Water'):
    initialize(substance)                                   #Initialize ports, files, and constants (see below) 
    temps = arange(startTemp, endTemp + incTemp, incTemp)   #Set temperatures to be used
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

def initialize(substance):#Method initializes motor, bath, file, and constants to be used in calculation
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
    line2 = 'Temp(C)\tShear Rate(1/s)\tShear Stress\tViscosity\tRotor Frequencies (1/s)\tMotor Frequency (1/s)'
    header = line1 + line2
    return header

def writeToFile(line):#Writes a string to file (called by main, initFile, closeUp)
    file = open(Filename, 'a')  #Open file to append
    file.write(line + '\n')     #Write line + escape char
    file.close()                #Close file so it can be accessed during data runs

def setConsts():#Sets constants for use in calculation (called by initialize)
    global GeoConst             #Initialize global GeoConst - geometrical constant for use in calculation
    global CalConst             #Initialize global CalConst - calibration constant for use in calculation
    ri = 11.0                   #Outer radius of the inner tube in mm
    ro = 12.95                  #Inner radius of the outer tube in mm
    GeoConst = 8.0 * math.pi * ro**2.0 * ri**2.0 * log(ro / ri) / (ro**2.0 - ri**2.0)**2.0 #Compute GeoConst
    CalConst = 2.0*5.965252117              #Set CalConst - must be recalibrated if tubes are switched

def setBath(temp, stirFreq):#Sets bath temperature and waits for equilibration (called by main)
    TBath.setSet(temp)          #Set setpoint temperature (see bath.py)
    q = setMotorFreq(stirFreq)  #Set a low frequency shear while temperature equilibrates
    print 'Waiting for temp to stabilize at ' + str(temp) + ' degrees C.'
    unstable = True
    minStableTime = 120.0       #Minimum time required for the bath to be at steady stayed before proceeding
    while unstable:
        t = time.time()
        unstable = False
        while abs(time.time() - t) < minStableTime: #Waits for stability for minStableTime
            if abs(TBath.extTemp() - temp) > .01:#Checks for temperature equal to the setpoint
                unstable = True
                break           #Restarts count if the bath temperature is not equal to setpoint at any point during the count
    print 'Temp = ' + str(temp) + ' degrees C. Stable setpoint reached.'#If minStableTime is reached, a message is printed and the program moves on

def motorCalibration(targMin, targMax, numMotorFreqs, motorCalFreq):#Calibrates the viscometer so that the minimum shear rate is close to targMin, the maximum shear rate is close to maxSR, and the others are equally spaced in between (called by main)
    rotorCalFreq = counter.counter(8)* math.pi / 180#Calibration Measurement (numerical factor 2pi/100 radians per spoke)
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

def processData(temp, motorFreq, rotorFreqs):#Computes data (visc, ss, sr) from angular frequencies of inner tube and angular frequency of motor, concatenates data into string, and returns data (called by main)
    motorFreqRadians = motorFreq * math.pi / 180.0          #Convert motor speed from degrees per second to radians per second
    rotorFreq = average(rotorFreqs)                         #Average all measurements
    v = CalConst *(motorFreqRadians - rotorFreq) / (rotorFreq)#Calculate viscosity
    ss = CalConst * GeoConst * (motorFreqRadians - rotorFreq)#Calculate shear stress
    sr = GeoConst * rotorFreq                               #Calculate shear rate 
    data = str(temp) + '\t' + str(sr) + '\t' + str(ss) + '\t' + str(v) + '\t' + str(rotorFreqs) + '\t' + str(motorFreq)#Create data String for writing to file
    return data                                             #Returns data string

def getData(numMeas):#Goes through all shear rates and samples each one numMeas times; returns angular frequencies of inner tube (called by main)
    rotorFreqs = zeros(numMeas)                             #Create array in which to store rotorFreqs
    for j in range (0, numMeas):                            #Take numMeas readings
        while rotorFreqs[j] == 0:                           #It is possible to have an unsuccessful count. This checks if any previous counts have been successful and allows another count if necessary
            rotorFreqs[j] = counter.counter(8)* math.pi / 180#Angular frequency measurement of inner tube (numerical factor 2pi/100 radians per spoke)
        print 'count ' + str(j + 1) + ' of ' + str(numMeas) + ' taken successfully.'
    return rotorFreqs                                       #Return angular frequencies of inner tube
            
def closeUp():#Ends Experiment (called by main)
    print 'Experiment Complete.'
    endTime = time.strftime("%y-%m-%d-%H-%M-%S", time.localtime())#System time
    writeToFile(endTime)                                    #Prints end time to file
    TBath.close()                                           #Closes temperature bath
    Motor.write('DI\r\n')                                   #Disables motor (see Manual.doc for details)
    Motor.close()                                           #Closes motor
    
def setMotorFreq(motorFreq):#Checks to see if motorFreq is small enough to be safe, sets the motor frequency, and waits for rotor to speed up (called by main, setBath)
    skipBool = False
    if motorFreq > 2000.0:                                  #Check motorFreq is safe
        print str(motorFreq) + ' deg/s too large. Measurement(s) skipped.'
        skipBool = True                                     #Set skipBool to True: tell program to skip data analysis for this measurement
    else:
        Motor.write('FR ' + str(motorFreq) + '\r\n')        #If motorFreq is safe, set motor to freerun at motorFreq (see Manual.doc for details on setting motor speed)
        print 'Motor frequency = ' + str(motorFreq) + ' degrees per second.'
        time.sleep(3)                                      #Wait for rotor to speed up
    return skipBool                                         #Returns skipBool that tells main to skip data analysis