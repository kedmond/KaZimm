# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 16:43:31 2013

@author: Kazem V. Edmond

Find inflection point in Zimm data.
Low shear rate data is linear, and Newtonian, but at some point we reach 
the Taylor Instability, resulting in a sudden change in slope.  Finding
that point automatically, and then fitting to only the Newtonian portion
would be nice.

The following routines find the inflection points in 2 different ways:
    1) Fitting a spline gives us the slope between each interval, and we 
    fit over the range where the slope does not vary much.
    2) Removing the mean trend of the data, we find the maximal peak of the
    transformed data, which corresponds to the inflection point.

This routine fits the data with a spline, using the built in function
within SciPy's interpolate library.  The spline gives us the slopes for
each increment.  The instability is marked by an abrupt change in slope.

"""

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline 

def findInflection1(x, y, threshold=0.9):
    """
    Fits y = m*x to Zimm data before onset of Taylor instability.
    Returns m.
    Last updated by Kazem Edmond on Feb. 26, 2013.
    
    Inflection is found by identifying sudden change in slope.  Slope
    between each interval is calculated using a spline fit.  Data is 
    cleaned by removing outliers using a simple "threshold", where
    difference in between adjacent x is within 5 times of first interval.
    The threshold parameter defines size of inflection to look for.
    
    """
    # Clean data by removing outliers to help spline fit
    wout = np.where( np.abs(x[0:-1] - x[1:]) > 5*np.abs(x[0]-x[1]) )
    wout2 = wout[0][range(1, np.size(wout), 2)]
    # Note: only use odd indices from np.where()
    
    # If necessary, new arrays for cleaned data:
    if np.size(wout2 > 0):
        xcln = np.delete(x, wout2)
        ycln = np.delete(y, wout2)
    else:
        xcln = x
        ycln = y
    
    # Fit a spline to cleaned data, giving slopes over intervals:
    spl = InterpolatedUnivariateSpline(xcln, ycln)
    
    # Store slopes in new array:
    nvals = int(np.floor(max(xcln)))
    spVals = np.zeros([nvals, 4])
    
    for i in range(1, nvals, 1): spVals[i] = spl.derivatives(i)
    
    # Find range over which it's Newtonian.
    # Only consider data that is within threshold of initial slope.
    wlin = np.where(spVals[1:, 1] > threshold * spVals[1, 1])
    
    wid = np.where(np.subtract(wlin, range(0, np.size(wlin)))[0] > 0)

    # Check if all of the data is Newtonian:    
    if np.size(wid) > 0:
        idd = int(wid[0])
        xlin = xcln[0:idd]
        ylin = ycln[0:idd]
    else:
        xlin = xcln[wlin]
        ylin = ycln[wlin]
    
    # 2.27.2013: Alternative to spline fitting
    # Remove overall trend of data, resulting peak is inflection point.
    # Requires that inflection peak actually exists
    # m, _, _, _ = np.linalg.lstsq(x[:, np.newaxis], y)
    # xp = x[yp == np.max(y - m * (x - np.mean(x)))]
    # mf, _, _, _ = np.linalg.lstsq(x[x<xp][:, np.newaxis], y[x<xp])
    
    # Least squares linear fit to get the slope:
    m, _, _, _ = np.linalg.lstsq(xlin[:, np.newaxis], ylin)
    
    return m

def findInflection2(x, y, threshold=0.9):
    """
    Fits y = m*x to Zimm data before onset of Taylor instability.
    Returns m.
    Last updated by Kazem Edmond on Feb. 26, 2013.
    
    Inflection is found by identifying sudden change in slope.  Slope
    between each interval is calculated using a spline fit.  Data is 
    cleaned by removing outliers using a simple "threshold", where 
    difference in between adjacent x is within 5 times of first interval.
    The threshold parameter defines size of inflection to look for.
    
    """
   
    # 2.27.2013: Alternative to spline fitting
    # Remove overall trend of data, resulting peak is inflection point.
    # Requires that inflection peak actually exists
    m, _, _, _ = np.linalg.lstsq(x[:, np.newaxis], y)
    yp = y - m * (x - np.mean(x))
    
    if np.std(yp) > 0.2:
        xp = x[yp == np.max(yp)]
        mf, _, _, _ = np.linalg.lstsq(x[x<xp][:, np.newaxis], y[x<xp])
        
        return mf
    else:
        return m

