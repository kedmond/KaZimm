'''
Created on Mar 9, 2010

This module reads the data file put out by Zimm.main and puts it into a 4-column matrix
with the temperature in the first column, viscosity in the second, shear stress in the 
third, and shear rate in the fourth:
temp    visc    ss    sr

This module is also programmed to plot temperature versus viscosity.

@author: bernie

'''
from pylab import *
import matplotlib.pyplot as plt
def an(filename):
    file = open(filename)        #opens file
    str = file.read()            #reads data out of file as giant string
    lines = str.split('\n')      #makes an array with each row of data as an entry, each row being a data point
    v = []                       #makes empty arrays for each type of information recorded
    ss = []                         
    sr = []
    t = []
    for line in lines[3:]:       #for array entries representing the fourth row onwards    #line = a row as an object
        nums = line.split('\t')  #each row (each entry in lines) split into individual pieces of info                  #nums = an array composed of each entry of the row as objects

        if len(nums) != 6: break #break if array nums doesn't have 6 entries, i.e. there arent 6 pieces of info

        nums.pop()               #pop returns the last item in nums and removes it
        nums.pop()               #these two pops discard rotor freq and motor freq     

        v.append(nums.pop())     #add next four popped nums into previously constructed arrays for each info piece type
        ss.append(nums.pop())    #data in file must be in order t, ss, sr, v, rf, mf (as it is in main's data output)
        sr.append(nums.pop())    
        t.append(nums.pop())
    data = [t, v, ss, sr]        #make a data array of the arrays of each piece of info
    plot(t, v, 'o')              #graph data (temp v viscosity) and label (default plot uses circle markers)
    xlabel('Temperature (C)')
    ylabel('Viscosity (proportional to Pa s)')
    title('Viscosity vs Temperature for SY FCF 27%')
    show()
    return data


    
    
