# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 15:22:52 2014

Saves Zimm viscometer data to an HDF5 file.  

Each data set will contain some number of columns of tick counts from the rotor
along with attributes describing the parameters used in the experiment.

Parameters:
    temperature
    depletant
    depletion concentration
    particle volume fraction
    background liquid
    liquid density
    NaCl concentration
    TMAH concentration
    F108 concentration
    rotor radius
    rotor height
    stator radius

To do:
    -create a HDF5 database from existing Zimm data directories


@author: kedmond
"""

import h5py as hdf
import glob

import numpy as np

def toHDF(directory):
    
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
    
    sampleAttr = np.chararray((paraNum, 2), itemsize=24)
    
    for i in range(paraNum):
        sampleAttr[i, 0] = parameterNames[i]
        sampleAttr[i, 1] = parameterList[i]
    
    #phInd = hdr[0].find('phi')
    #phiStr = hdr[0][phInd+3:phInd+9].strip()
    # name of the HDF file we are creating...maybe use an existing file?
    #dateDirs = glob.glob(directory+'*2014')         # all the experiments were done this year
    
    #for dir in dateDirs:
    #    dateStr = dir[dir.find('/', len(dir)-14)+1:]  # dates are less than 12 chars long
    #    filename = 'database-'+dateStr+'.hdf5'
    
    # let's assume that we're looking at a more current directory:
    dateStr = directory[directory.find('/', len(directory)-14)+1:]
    filename = 'database-'+dateStr+'.hdf5'
    
    # create a file, for now
    outFile = hdf.File(filename)
    
    # Load in the data
    # Experiment directories
    expDirs = glob.glob(directory+'/Experiment*')
    expDirs.sort()
    
    # the experiments need attributes too
    
    # For each experiment:
    #   load and parse the descriptionfile
    #   load and parse each of the rotorSpeed files
    for exp in expDirs:
        # load up the file's data and experimental attributes
        descriptor = glob.glob(exp+'/description*')[0]
        descriptData = np.genfromtxt(descriptor, delimiter=',')
        
        dset = outFile.create_dataset(dateStr, data=descriptData, dtype='i')
        
        # get parameters from the file's header
        f = open(descriptor)                     # open the description file
        hdr = [f.readline(), f.readline()] 	 # first two lines contain info
        f.close()
        
        rotorFiles = glob.glob(exp+'/rotor*')
        rotorFiles.sort()
        
        for inFile in rotorFiles:
            # the first row is the calculated frequency from the FFT
            data = np.genfromtxt(inFile, delimiter=',')
            
        
    
    # create a dataset, with our newly loaded data
    outFile.create_dataset()
    
    outFile.close()
    
    
    