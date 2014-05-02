# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 23:40:58 2014

This routine digs down into a given directory, finds data sets labeled "rotorSpeeds*"
and plots them as eta/C vs. gamma dot (relative viscosity vs. shear rate)

A lot of things are assumed here, but it works well.

The directory can be anywhere, but the files must be named "rotorSpeeds*" and
have the data stored as comma-separated columns, where the data are the counts
from the Zimm's rotor.

To do items:
    -I need it to plot hysterisis as well.  That's a nice feature.
    -It would be nice to plot multiple "Experiments" in a single plot, perhaps
    by accepting multiple directories.

@author: kedmond
"""

import matplotlib.pyplot as plt
import brewer2mpl
import KaZimm5 as KZ
import numpy as np
import glob


def plotData(directory):
    
    # load in the data
    wmwr = KZ.loadData(directory)
    
    # get the phi
    with open(glob.glob(directory+'description*')[0], 'r') as f: sampleStr = f.readline()
    phiStr = sampleStr[sampleStr.find('phi')+3:sampleStr.find('phi')+9].strip()
    
    # Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
    set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
    
    # Save a nice dark grey as a variable
    almost_black = '#262626'
    
    fig, ax = plt.subplots(1)
    fig.patch.set_facecolor('white')
    fig.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.15)
    plt.grid()
    
    plt.ylabel(r'$\frac{\eta}{C}$', fontsize=28)
    plt.xlabel(r'$\dot{\gamma} \, (s^{-1})$', fontsize=24)
    
    # Plot the data here
    
    # shear rate and relative viscosity/calibration constant
    data = wmwr[:][wmwr['wr']>0]
    
    x = 12 * data['wr']
    y = (data['wm'] - data['wr'])/data['wr']
    
    
    ax.plot(x, y, color=set2[0], linewidth=0.75, label=r'$\phi = $'+phiStr)
    ax.semilogx()
    ax.semilogy()
    
    # plot each point with a unique color, corresponding with motor speed
    for i in range(x.size):
        color = set2[i]
        ax.scatter(x[i], y[i], label=str(int(data['wm'][i])), edgecolor='none', facecolor=color, s=100)
        
    # Error bars, from wmwr['err']
    
    
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
    
    
    # The following is for drawing a nice legend, which I am not using just yet
    
    # Remove the line around the legend box, and instead fill it with a light grey
    # Also only use one point for the scatterplot legend because the user will 
    # get the idea after just one, they don't need three.
    light_grey = np.array([float(248)/float(255)]*3)
    legend = ax.legend(frameon=True, scatterpoints=1)
    rect = legend.get_frame()
    rect.set_facecolor(light_grey)
    rect.set_linewidth(0.0)
    
    # Change the legend label colors to almost black, too
    texts = legend.texts
    
    for t in texts:
        t.set_color(almost_black)
    
    #ax.set_title(directory)
    #fig.savefig('scatter_matplotlib_improved_10_pretty_legend.png')


