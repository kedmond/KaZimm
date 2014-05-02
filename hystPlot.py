# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 13:07:19 2014

@author: kedmond
"""


import matplotlib.pyplot as plt
import numpy as np
import brewer2mpl

def plot(directory, filename='hystPlot.pdf'):

    # Load data from CSV to Numpy array:
    data = np.genfromtxt(directory, skip_header=2, delimiter=',', usecols = (1, 2), dtype=('float64, float64'))
    data.dtype.names = ('wm', 'wr')    # label the columns
    
    data['wm'] *= 180/np.pi
    
    x = 12*data['wr']
    y = (data['wm'] - data['wr']) / data['wr']
    
    # Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
#    set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
        
    # Save a nice dark grey as a variable
    almost_black = '#262626'
    
    fig, ax = plt.subplots(1)
    fig.patch.set_facecolor('white')
   #fig.set_size_inches(1025/72.0, 400/72, forward=True)
    
    # Set the axes titles
    #ax.set_title(str(motorFreq) + ' deg/sec', fontsize=24)
    
    # set up axes labels with greek symbols:
    ax.set_xlabel(r'$\dot{\gamma} \, (s^{-1})$', fontsize=24)
    ax.set_ylabel(r'$\frac{\eta}{C}$', fontsize=28)

    #fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.15)
    
    ax.grid(True)    
    
    # Plot the data here
    ax.plot(x, y, color='black')
    ax.semilogx()
    ax.semilogy()
    
    ax.plot([10, 1e4+1], [0.7, 0.7], '--', linewidth=4, color='g')
    
    ax.set_ylim([0.5, 100])
    ax.set_xlim([10, 1e4+100])
    
    # 11 RGB color values from ColorBrewer
    colorSet = brewer2mpl.get_map('Paired', 'qualitative', 11).mpl_colors
    
    
    
    # Remove top and right axes lines ("spines")
    spines_to_remove = ['top', 'right']
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)
    
    # Get rid of ticks. The position of the numbers is informative enough of
    # the position of the value.
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    
    # For remaining spines, thin out their line and change the black to a slightly off-black dark grey
    almost_black = '#262626'
    spines_to_keep = ['bottom', 'left']
    for spine in spines_to_keep:
        ax.spines[spine].set_linewidth(0.5)
        ax.spines[spine].set_color(almost_black)
    
    # Change the labels to the off-black
    ax.xaxis.label.set_color(almost_black)
    ax.yaxis.label.set_color(almost_black)
    
    # Change the axis title to off-black
    ax.title.set_color(almost_black)
    
    if filename != '':
        fig.savefig('/home/kedmond/Dropbox/ZimmViscometer/Plots/BTF2 - phi0.133 cp0.7/hysterisisPlots/'+
        filename, format='pdf')

    return (fig, ax)
