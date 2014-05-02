# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:22:55 2012

@author: kedmond

Written on 10/11/2012

To Do: Want to make this a subroutine for KaZimm
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

# file location
d = 'E:\Dropbox\Zimm viscometer'
f = '12-09-28-15-50-39_NaBr+H2O.txt'


def plotData(dir=d, filename=f):  # plot the data with coefficients

     # load up the data
    if sys.platform == 'linux2':
        data = np.genfromtxt(dir+'/'+filename, delimiter=',', skip_header=2, usecols=(0,1,2))
    if sys.platform == 'win32':
        data = np.genfromtxt(dir+'\\'+filename, delimiter=',', skip_header=2, usecols=(0,1,2))
    #stderr = np.genfromtxt(dir+filename, delimiter=',', skip_header=2, usecols=(3))
    
    # Get number of temps for plot range:
    tval = len(np.where(data[:, 0] == data[0, 0])[0])
    # The above assumes same amount of data for each temperature.
    
    # How many temperatures?
    temps = np.unique(data[:, 0])
    tnum = len(temps)
    
    # Store the slopes here
    clist = np.zeros(tnum, dtype='float32')
    
    # White background >> gray background
    plt.figure(facecolor='w')
    
    plt.plot()
    
    # set up axes labels
    plt.xlabel(r'$\omega_m - \omega_r$ ($\frac{\rm{rad.}}{s}$)', fontsize=20)
    plt.ylabel(r'$\omega_r$ ($\frac{\rm{rad.}}{s}$)', fontsize=20)
    
    for i in range(0, tnum, 1):
        # temp variable for the data
        x = data[i*tval:(i+1)*tval, 1] - data[i*tval:(i+1)*tval, 2]
        y = data[i*tval:(i+1)*tval, 2]
    
        # calculate the coefficient
        clist[i], _, _, _ = np.linalg.lstsq(x[:, np.newaxis], y)
        
        # plot data
        plt.plot(x, y, 'o', label='T = %d$^\circ$, $C$ = %f' % (temps[i], clist[i]))
        
        # plot the line of best fit
        plt.plot(x, clist[i]*x, '-r')
        
    
    # draw a legend:
    plt.legend(loc='lower right')
