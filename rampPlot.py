# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 09:08:11 2014

@author: kedmond
"""

import matplotlib.pyplot as plt
import numpy as np

def plot(directory, color, xlims, filename='rampPlot.pdf'):

    # get the parameters
    f = open(directory, 'r')
    hdr = f.readline()
    f.close()
    
    parameters = [int(s) for s in hdr.replace(',', '').split() if s.isdigit()]
    
    #duration = parameters[0]
    sampRate = parameters[1]
    #motorFreq = parameters[2]
    
    data = np.genfromtxt(directory, skiprows=2, delimiter=', ')
    
    # Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
#    set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
    
    # Save a nice dark grey as a variable
    almost_black = '#262626'
    
    fig, ax = plt.subplots(1)
    fig.patch.set_facecolor('white')
    fig.set_size_inches(1025/72.0, 400/72, forward=True)
    
    # Set the axes titles
    #ax.set_title(str(motorFreq) + ' deg/sec', fontsize=24)
    
    ax.set_xlabel('duration (s)', fontsize=24)
    
    fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.15)
    
    ax.grid(True)    
    
    # Plot the data here
    ax.plot(data[:, 0], linewidth=0.75, color=color)
    ax.set_xlim(xlims)
    ax.set_ylim([0, 11])
  
    # 28.Feb.2014, Kaz: change x-axis to seconds
    xlabels = ax.get_xticks()
    xlabels /= sampRate
    ax.set_xticklabels(np.round(xlabels, 1), fontsize=24)
    
    ax.set_yticklabels([])

  
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
        fig.savefig('/home/kedmond/Dropbox/ZimmViscometer/Plots/BTF2 - phi0.133 cp0.7/rampPlots/'+
        filename, format='pdf')

    return (fig, ax)
