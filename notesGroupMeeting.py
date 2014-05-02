# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 11:39:03 2014

Notes for today's group meeting

@author: kedmond
"""

import rampPlot as rp
import hystPlot as hp
import brewer2mpl

directory = '/media/HDD/Data/ZimmViscometer/Data/BTF2/17Apr2014/Experiment9/'

set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors

# plots for ramping up
rp.plot(directory+'rotorSpeedsUp-0200.txt', set2[0], [0, 1000], 'rampUp-0200.pdf')
rp.plot(directory+'rotorSpeedsUp-0600.txt', set2[1], [0, 2000], 'rampUp-0600.pdf')
rp.plot(directory+'rotorSpeedsUp-1000.txt', set2[2], [0, 6000], 'rampUp-1000.pdf')
rp.plot(directory+'rotorSpeedsUp-2000.txt', set2[3], [0, 10000], 'rampUp-2000.pdf')

# plots for ramping down
rp.plot(directory+'rotorSpeedsDown-0200.txt', set2[0], [0, 1000], 'rampDown-0200.pdf')
rp.plot(directory+'rotorSpeedsDown-0600.txt', set2[1], [0, 2000], 'rampDown-0600.pdf')
rp.plot(directory+'rotorSpeedsDown-1000.txt', set2[2], [0, 6000], 'rampDown-1000.pdf')
rp.plot(directory+'rotorSpeedsDown-2000.txt', set2[3], [0, 10000], 'rampDown-2000.pdf')

# plot hysterisis

directory = '/media/HDD/Data/ZimmViscometer/Data/BTF2/9.April.2014/Set 1/'
filename = '14-04-09-09-32-44_BTF2 - phi 0.133 - 0.7 g_L PEO 600k, D2O_H2O, 1.053 g_mL.txt'

hp.plot(directory+filename, [0, 1000], 'hysterisis.pdf')
