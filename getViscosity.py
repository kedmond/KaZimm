# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 14:18:38 2013

@author: Kazem V. Edmond

This program will calculate the viscosity of a fluid based on Zimm-Crothers
viscometer data.

History:
    Updated by KVE, 3.06.2013
"""

import numpy as np



def getViscosity(wr_0, wr, C, eta_0 = 0.0014):
    # viscosity of NaBr in H2O is 1.4 mPa*s at 20C (tabulated)
    # viscosity of H2O is 1 mPa s at 20C, 0.89 mPa s at 25C, 0.797 mPa s at 30C
    wr_0mean = np.mean(wr_0)
    wrmean = np.mean(wr)
    w = wr_0mean/wrmean
    
    viscosity = w*eta_0 + C*(w-1)
    
    return viscosity

