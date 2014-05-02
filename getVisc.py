# -*- coding: utf-8 -*-
"""
Created on Sat Feb 23 13:47:04 2013

@author: kedmond

Routine for calculating viscosity from wr, wm, C, and eta_0.

User gives it a row from a new data set, [T, wm, wr], and the full data 
set from the original background liquid.  The routine will match up T and
wm to find wr0, and calculate the sample's viscosity from that.

"""

import numpy as np
import matplotlib.pyplot as plt
import os                       # detect OS being used
import subprocess as subproc    # used to open PDF in Linux
import time                     # create timestamps in filenames


# sampledata = [T, wm, wr]

def getVisc(bgvisc, bgdata, sampledata):
    
    # insert math here
    
    return viscosity
    




