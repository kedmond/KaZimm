'''
Created on Jun 25, 2009

@author: bernie
''' 
import serial

class RW(object):
    def __init__(self, port='COM5', to=0.1):
        self.port = port
        self.bathPort = serial.Serial(port, timeout=to)
    
    #write to bath    
    def write(self, str):
        self.bathPort.write(str)
    
    #Read from bath    
    def read(self):
        dataList = []
        data = self.bathPort.read()
        while data != '':
            dataList = dataList + [ord(data)]
            data = self.bathPort.read()
        return dataList
    
    #Retrieve Internal Temp in degrees C    
    def intTemp(self):
        str = '' + chr(0xCA) + chr(0x00) + chr(0x01) + chr(0x20) + chr(0x00) + chr(0xDE)
        self.write(str)
        data = self.read()
        temp = 25.6 * data[6] + 0.1*data[7]
        if temp > 200:
            temp = temp - 6553.6
        return temp
    
    #Retrieve External Sensor Temp in degrees C
    def extTemp(self):
        str = '' + chr(0xCA) + chr(0x00) + chr(0x01) + chr(0x21) + chr(0x00) + chr(0xDD)
        self.write(str)
        data = self.read()
        temp = 25.6 * data[6] + 0.1*data[7]
        if temp > 200:
            temp = temp - 6553.6
        return temp
    
    #Retrieve Setpoint Temp in degrees C
    def getSet(self): 
        str = '' + chr(0xCA) + chr(0x00) + chr(0x01) + chr(0x70) + chr(0x00) + chr(0x8E)
        self.write(str)
        data = self.read()
        temp = 25.6 * data[6] + 0.1*data[7]
        if temp > 200:
            temp = temp - 6553.6
        return temp
    
    #Set Setpoint Temp in degrees C
    def setSet(self, stpt):
        if stpt<0:
            stpt = stpt + 6553.6
        stpt = 10*stpt
        d1 = int(stpt/256)
        d2 = int(stpt%256)
        sum = 0x00 + 0x01 + 0xF0 + 0x02 + d1 + d2
        cs = int(255 ^ (sum%256))
        str = '' + chr(0xCA) + chr(0x00) + chr(0x01) + chr(0xF0) + chr(2) + chr(d1) + chr(d2) + chr(cs)
        self.write(str)
        data = self.read()
        temp = 25.6 * data[6]# + 0.1*data[7]
        if temp > 200:
            temp = temp - 6553.6
        return temp
        
    #close the bath
    def close(self):
        self.bathPort.close()
    
            
