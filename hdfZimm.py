# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 15:22:52 2014

Saves Zimm viscometer data to an HDF5 file.  

Each data set will contain some number of columns of tick counts from the rotor
along with attributes describing the parameters used in the experiment.

Parameters:
    temperature
    depletion type
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

@author: kedmond
"""

import h5py as hdf
import numpy as np



def writeData(data):
    
    f = hdf.File('mytestfile.hdf5', 'w')
    
    dset = f.create_dataset('mydataset', (100,), dtype='i')
    
    dset[...] = np.arange(100)
    
    dset[0]
    dset[10]
    dset[0:100:10]
    
    dset.name
    
    f.name
    
    grp = f.create_group('subgroup')
    
    dset2 = grp.create_dataset('another_dataset', (50,), dtype='f')
    
    dset2.name
    
    dset3 = f.create_dataset('subgroup2/dataset_three', (10,), dtype='i')
    dset3.name

    f['subgroup2/dataset_three']

    for name in f:
        print name
        
    dset.attrs['temperature'] = 99.5
    dset.attrs['temperature']
    