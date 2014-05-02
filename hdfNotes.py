# -*- coding: utf-8 -*-
"""
Created on Thu May  1 16:22:02 2014

@author: kedmond
"""

import h5py as hdf
import glob

import numpy as np


directory = '/media/HDD/Data/ZimmViscometer/Data/BTF2/'
dateDirs = glob.glob(directory+'*Apr2014')         # all the experiments were done this year
dateDirs.sort()
dateStr = dateDirs[2]
dateStr = dateStr[dateStr.find('/', len(dateStr)-14)+1:]

filename = 'database-'+dateStr+'.hdf5'
outFile = hdf.File(filename)

expDirs = glob.glob(directory + '/' + dateStr + '/Experiment*')
expDirs.sort()
expDir = expDirs[2]
expStr = expDir[expDir.find('/', len(expDir)-14)+1:]

dataFiles = glob.glob(expDir + '/rotor*txt')
dataFiles.sort()
dataName = dataFiles[0][dataFiles[0].find('/', len(dataFiles[0])-30)+1:]

data = np.genfromtxt(dataFiles[0], delimiter=',', skiprows=2)

grp = outFile.create_group(expStr)
dset0050 = grp.create_dataset(dataName, data=data)

# Prepare the list of attributes (experimental parameters)
parameterNames = ['date', 'experiment type', 'temperature', 'depletant', 
                 'depletant concentration', 'particle volume fraction', 
                 'liquid', 'liquid density', 'NaCl concentration', 
                 'TMAH concentration', 'F108 concentration', 'rotor radius', 
                 'rotor height', 'stator radius']
parameterList = ['empty', 'empty', '20.0', 'PEO 600k', 'empty', 'empty',
                 'D2O_H2O', '1.053 g/mL', '10 mM', '0.01wt%', '0.05wt%', 
                 '12 mm', '3 inches', '13 mm']
                 
paraNum = np.size(parameterNames)

for i in range(paraNum):
    grp.attrs[parameterNames[i]] = parameterList[i]

dataParameterNames = ['motor frequency', 'sample rate', 'duration', 'motor time', 'pause time']
dataParameters = ['50 deg/s', '16 counts/s', '600 s', '60 s', '120 s']

for i in range(np.size(dataParameterNames)):
    dset0050.attrs[dataParameterNames[i]] = dataParameters[i]

