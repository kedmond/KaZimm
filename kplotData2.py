# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:22:55 2012

@author: Kazem V. Edmond

Written on 10/11/2012
Updated on 12/19/2012 by KVE: sorts data by temperature.
Updated on 2/26/2013 by KVE: fits before inflection point.

To Do: Want to make this a subroutine for KaZimm
"""

import numpy as np
import matplotlib.pyplot as plt
import os                       # detect OS being used
import subprocess as subproc    # used to open PDF in Linux
import time                     # create timestamps in filenames

# Use spline fitting to find Newtonian region:
from scipy.interpolate import InterpolatedUnivariateSpline 

# plot the data with coefficients
def plotData(directory, filename, fit=True, writeFile=False):  
    """
    Function for plotting data from the Zimm viscometer.
    Last updated by Kazem Edmond on Dec. 19, 2012.    
    
    The Zimm measures the angular speed of the magnets omega_m and of the 
    rotor in the sample omega_r.    
    This function plots omega_r versus the difference, omega_m - omega_r.
    The slope C is the calibration constant for that particular sample at
    the given temperature, which can be used to determine the relative
    viscosity eta_r.
    """
    if os.name == 'posix':
        slstr = '/'
    elif os.name == 'nt':
        slstr = '\\'
    
    # Load data from CSV to Numpy array:
    data = np.genfromtxt(directory+slstr+filename, skip_header=2, delimiter=',', 
                         usecols = (0, 1, 2), dtype=('float64, float64, float64'))
    data.dtype.names = ('T', 'wm', 'wr')    # label the columns
    
    # sort in terms of temperature:
    data = np.sort(data, order='T')
    
# Uncomment for errorbars...unlikely.
#    # some data has 4th column with std. error:
#    if len(data[0, :]) == 4:
#        stderr = data[:, 3:4]
#        data = data[:, 0:3]
#    elif len(data[0, :]) == 3:
#        data = data[:, 0:3]
#        stderr = np.zeros([1])

    # Count and store temperatures (assume equal number of each temp):
    tval = len( np.where(data['T'] == data['T'][0])[0] )
    temps = np.unique(data['T'])
    tnum = len(temps)
    
    if fit:
        clist = np.zeros(tnum, dtype='float32')     # Store the slopes here

    # Setup a figure with a new plot:
    plt.figure(facecolor='w') # Gray doesn't work for pubs or talks.
    plt.hold(True)    # Hold figure open to plot repeatedly to same axes
    plt.clf()         # Clear the plot before adding new data
   
    # set up axes labels with greek symbols:
    plt.xlabel(r'$\omega_m - \omega_r$ ($\frac{\rm{rad.}}{s}$)', fontsize=20)
    plt.ylabel(r'$\omega_r$ ($\frac{\rm{rad.}}{s}$)', fontsize=20)
    
    # Loop through data from each Temp., plotting to same axes:
    for i in range(0, tnum, 1):
        # temp variable for the data
        x = data['wm'][i*tval:(i+1)*tval] - data['wr'][i*tval:(i+1)*tval]
        y = data['wr'][i*tval:(i+1)*tval]

        if fit:
            # calculate the coefficient
            # clist[i], _, _, _ = np.linalg.lstsq(x[:, np.newaxis], y)
            clist[i] = findInflection(x, y)
            

            # Plot data for this temperature, w/ fit:
            plt.plot(x, y, 'o', label='T = %d$^\circ$, $C$ = %f' % (temps[i], clist[i]))
            # plot the line of best fit
            plt.plot(x, clist[i]*x, '-r')
        else:
            # Plot data for this temperature w/o fit:
            plt.plot(x, y, 'o', label='T = %d$^\circ$' % (temps[i]))


# Uncomment if you want errorbars...unlikely:
#        if stderr.any():
#            yerr = (stderr.flatten())[i*tval:(i+1)*tval]
#            plt.errorbar(x, y, yerr, fmt='o', label='T = %d$^\circ$, $C$ = %f' % (temps[i], clist[i]))
#        else:
#            plt.plot(x, y, 'o', label='T = %d$^\circ$, $C$ = %f' % (temps[i], clist[i]))
    
    # draw a legend:    
    plt.legend(loc='lower right')
    plt.grid(True)
    
    plt.hold(False)
    
    # save the plot in a folder, remove .txt from filename:
    if writeFile:
        # create unique name of data directory:
        today = time.strftime("%y-%m-%d", time.localtime())
        plotsPath = directory + slstr + 'plotDataKZ-plots'+slstr+today
        
        # string of file path:
        filePath = plotsPath + slstr + filename.replace(filename[-3:], 'pdf')
       
        # if need be, make a new directory:
        if not os.path.exists(plotsPath):
            os.makedirs(plotsPath)
        
        # save the plot to data directory and open it with PDF reader:
        plt.savefig(filePath, format = 'pdf')
    
        if os.name == 'posix':
            subproc.call(('xdg-open', filePath))
        elif os.name == 'nt':
            os.startfile(filePath)

    else:
        # we're not creating a PDF, so just show in Python:
        plt.show()

    # It's important to close out variables:
    # See http://matplotlib.org/users/pyplot_tutorial.html for details.
    #plt.close()


def findInflection(x, y, threshold=0.9):
    """
    Fits y = m*x to Zimm data before onset of Taylor instability.
    Returns m.
    Last updated by Kazem Edmond on Feb. 26, 2013.
    
    Inflection is found by identifying sudden change in slope.  Slope
    between each interval is calculated using a spline fit.  Data is 
    cleaned by removing outliers using a simple "threshold", where difference
    in between adjacent x is within 5 times of first interval.
    The threshold parameter defines size of inflection to look for.
    
    """
    # Clean data by removing outliers to help spline fit
    wout = np.where( np.abs(x[0:-1] - x[1:]) > 5*np.abs(x[0]-x[1]) )
    wout2 = wout[0][range(1, np.size(wout), 2)]    # Note: only use odd indices from np.where()
    
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

