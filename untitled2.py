# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 15:18:27 2014

@author: kedmond
"""

import h5py as hdf
import numpy as np

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

dataset_three = f['subgroup2/dataset_three']

for name in f:
    print name


