# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:22:55 2012

@author: Kazem V. Edmond

Written on 10/11/2012
Updated on 12/19/2012 by KVE: sorts data by temperature.
Updated on 2/26/2013 by KVE: fits before inflection point.
Updated on 3/07/2013 by KVE: added 'xmax' parameter

Parameters:
    
    directory: String with directory path containing data
    filename: String with file name of data set
    writeFile[FALSE]: Boolean to write plot to PDF
    xmax[-1]: If defined, data greater than xmax will be cropped out

To Do:  Want to make this a subroutine for KaZimm
        March '13: Parameter to pick b/w spline() or trend fitting

"""

import numpy as np
import matplotlib.pyplot as plt
import os                       # detect OS being used
import subprocess as subproc    # used to open PDF in Linux
import time                     # create timestamps in filenames

# plot the data with coefficients
def plotData(directory, filename, fit=True, writeFile=False, xmax=-1):  
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

        # crop out data if user specifies
        # if cropped, that's what is fitted to below        
        if xmax > -1:
            x = x[x < xmax]
            y = y[x < xmax]
        
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

