# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:22:55 2012

@author: Kazem V. Edmond

Written on 10/11/2012

To Do: Want to make this a subroutine for KaZimm
"""
import numpy as np
import matplotlib.pyplot as plt
import os                       # detect OS being used
import subprocess as subproc    # used to open PDF in Linux
import time                     # create timestamps in filenames


# plot the data with coefficients
def plotData(directory, filename, writeFile=False):  
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
    
    # Load up the data:
    data = np.genfromtxt(directory+slstr+filename, delimiter=',', skip_header=2, usecols=(0,1,2))
    #data = np.genfromtxt(directory+slstr+filename, delimiter=',', skip_header=2)
    #data = np.genfromtxt(files[24], skip_header=2, delimiter=',', usecols = (0, 1, 2), dtype=('float64, float64, float64'))
    #data.dtype.names = ('Temps', 'Motor Speed', 'Rotor Speed')
    
    
# Uncomment for errorbars...unlikely.
#    # some data has 4th column with std. error:
#    if len(data[0, :]) == 4:
#        stderr = data[:, 3:4]
#        data = data[:, 0:3]
#    elif len(data[0, :]) == 3:
#        data = data[:, 0:3]
#        stderr = np.zeros([1])

    # Count and store temperatures:
    tval = len(np.where(data[:, 0] == data[0, 0])[0])
    # The above assumes same amount of data for each temperature.
    temps = np.unique(data[:, 0])
    tnum = len(temps)
    
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
        x = data[i*tval:(i+1)*tval, 1]-data[i*tval:(i+1)*tval, 2]
        y = data[i*tval:(i+1)*tval, 2]
        # calculate the coefficient
        clist[i], _, _, _ = np.linalg.lstsq(x[:, np.newaxis], y)
        
        # Plot data for this temperature, w/ or w/o error bars:
        plt.plot(x, y, 'o', label='T = %d$^\circ$, $C$ = %f' % (temps[i], clist[i]))

# Uncomment if you want errorbars...unlikely:
#        if stderr.any():
#            yerr = (stderr.flatten())[i*tval:(i+1)*tval]
#            plt.errorbar(x, y, yerr, fmt='o', label='T = %d$^\circ$, $C$ = %f' % (temps[i], clist[i]))
#        else:
#            plt.plot(x, y, 'o', label='T = %d$^\circ$, $C$ = %f' % (temps[i], clist[i]))
        
        
        # plot the line of best fit
        plt.plot(x, clist[i]*x, '-r')
    
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
        filePath = plotsPath + slstr + filename.replace('.txt', '', 1)+'.pdf'
       
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
    plt.close()
