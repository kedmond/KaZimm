# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 21:20:20 2014

@author: kedmond
"""

# plot of wr vs. wm for BTF2 particles in Zimm

import matplotlib.pyplot as plt
import numpy as np
import glob

dirname = '/home/kedmond/Dropbox/ZimmViscometer/Data/BTF2/'

files = glob.glob(dirname+'*.txt')

data1 = np.genfromtxt(files[2], delimiter = ',', skiprows=2)
data2 = np.genfromtxt(files[3], delimiter = ',', skiprows=2)
data3 = np.genfromtxt(files[4], delimiter = ',', skiprows=2)

wm = np.array(data1[:, 1]) * 180/np.pi

wr1 = np.array(data1[:, 2]) * 180/np.pi
wr2 = np.array(data2[:, 2]) * 180/np.pi
wr3 = np.array(data3[:, 2]) * 180/np.pi


fig = plt.figure(facecolor='white')

plt.loglog(wm, wr1, '-o')
plt.hold(True)
plt.loglog(np.append(wm, 1000), wr2, '-o')
plt.loglog(wm, wr3, '-o')

ranges = np.array([1, 1000])

plt.loglog(ranges, ranges*1.5, '--')
plt.loglog(ranges, ranges/2, '--')

plt.xlabel(r'$\omega_m (deg./sec.)$', fontsize=20)

plt.ylabel(r'$\omega_r (deg./sec.)$', fontsize=20)

plt.title('Zimm: 900 nm particles, 0.8 g/L PEO, 1.051 g/mL')

plt.grid(True)

plt.hold(False)


fig = plt.figure(facecolor='white')

plt.loglog(wm-wr1, wr1, '-o')
plt.hold(True)
plt.loglog(np.append(wm, 1000)-wr2, wr2, '-o')
plt.loglog(wm-wr3, wr3, '-o')

ranges = np.array([1, 1000])

plt.loglog(ranges, ranges*1.5, '--')

plt.xlabel(r'$\omega_m - \omega_r (deg./sec.)$', fontsize=20)

plt.ylabel(r'$\omega_r (deg./sec.)$', fontsize=20)

plt.title('Zimm: 900 nm particles, 0.8 g/L PEO, 1.051 g/mL')

plt.grid(True)

plt.hold(False)


# plot calibration curve

dirname = '/home/kedmond/Dropbox/ZimmViscometer/Data/Calibration/'

files = glob.glob(dirname+'*0.8*.txt')

calib1 = np.genfromtxt(files[0], delimiter = ',', skiprows=2)
calib2 = np.genfromtxt(files[1], delimiter = ',', skiprows=2)
calib3 = np.genfromtxt(files[2], delimiter = ',', skiprows=2)

wm = np.array(calib1[:, 1]) * 180/np.pi

wr1 = np.array(calib1[:, 2]) * 180/np.pi
wr2 = np.array(calib2[:, 2]) * 180/np.pi
wr3 = np.array(calib3[:, 2]) * 180/np.pi


fig = plt.figure(facecolor='white')

plt.loglog(wm-wr1, wr1, '-o')
plt.hold(True)
plt.loglog(wm-wr2, wr2, '-o')
plt.loglog(wm-wr3, wr3, '-o')

ranges = np.array([1, 1000])

plt.loglog(ranges, ranges*1.5, '--')

plt.xlabel(r'$\omega_m - \omega_r (deg./sec.)$', fontsize=20)

plt.ylabel(r'$\omega_r (deg./sec.)$', fontsize=20)

plt.title('Calibration: 0.8 g/L PEO, 1.051 g/mL')

plt.grid(True)

plt.hold(False)


fig = plt.figure(facecolor='white')

plt.loglog(wm, wr1, '-o')
plt.hold(True)
plt.loglog(wm, wr2, '-o')
plt.loglog(wm, wr3, '-o')

ranges = np.array([1, 1000])

plt.loglog(ranges, ranges*1.5, '--')

plt.xlabel(r'$\omega_m (deg./sec.)$', fontsize=20)

plt.ylabel(r'$\omega_r (deg./sec.)$', fontsize=20)

plt.title('Calibration: 0.8 g/L PEO, 1.051 g/mL')

plt.grid(True)

plt.hold(False)


# plot the calibration data with the chain data:

datadir = '/home/kedmond/Dropbox/ZimmViscometer/Data/BTF2/'

files = glob.glob(datadir+'*.txt')

data1 = np.genfromtxt(files[1], delimiter = ',', skiprows=2)
data2 = np.genfromtxt(files[2], delimiter = ',', skiprows=2)
data3 = np.genfromtxt(files[3], delimiter = ',', skiprows=2)

wm = np.array(data1[:, 1]) * 180/np.pi

wr1 = np.array(data1[:, 2]) * 180/np.pi
wr2 = np.array(data2[:, 2]) * 180/np.pi
wr3 = np.array(data3[:, 2]) * 180/np.pi


calibdir = '/home/kedmond/Dropbox/ZimmViscometer/Data/Calibration/'

files = glob.glob(calibdir+'*0.8*.txt')

calib1 = np.genfromtxt(files[0], delimiter = ',', skiprows=2)
calib2 = np.genfromtxt(files[1], delimiter = ',', skiprows=2)
calib3 = np.genfromtxt(files[2], delimiter = ',', skiprows=2)

wmcal = np.array(calib1[:, 1]) * 180/np.pi

wr1cal = np.array(calib1[:, 2]) * 180/np.pi
wr2cal = np.array(calib2[:, 2]) * 180/np.pi
wr3cal = np.array(calib3[:, 2]) * 180/np.pi


fig = plt.figure(facecolor='white')

plt.loglog(wm, wr1, '-o')
plt.hold(True)
plt.loglog(np.append(wm, 1000), wr2, '-o')
plt.loglog(wm, wr3, '-o')

plt.loglog(wmcal, wr1cal, '-o')
plt.loglog(wmcal, wr2cal, '-o')
plt.loglog(wmcal, wr3cal, '-o')


ranges = np.array([1, 1000])

plt.loglog(ranges, ranges*1.5, '--')
plt.loglog(ranges, ranges/2, '--')

plt.xlabel(r'$\omega_m (deg./sec.)$', fontsize=20)

plt.ylabel(r'$\omega_r (deg./sec.)$', fontsize=20)

plt.title('0.8 g/L PEO, 1.051 g/mL')

plt.grid(True)

plt.hold(False)


fig = plt.figure(facecolor='white')

plt.loglog(wm-wr1, wr1, '-o')
plt.hold(True)
plt.loglog(np.append(wm, 1000)-wr2, wr2, '-o')
plt.loglog(wm-wr3, wr3, '-o')

plt.loglog(wmcal-wr1cal, wr1cal, '-o')
plt.loglog(wmcal-wr2cal, wr2cal, '-o')
plt.loglog(wmcal-wr3cal, wr3cal, '-o')


ranges = np.array([1, 1000])

plt.loglog(ranges, ranges*1.5, '--')

plt.xlabel(r'$\omega_m-\omega_r (deg./sec.)$', fontsize=20)

plt.ylabel(r'$\omega_r (deg./sec.)$', fontsize=20)

plt.title('0.8 g/L PEO, 1.051 g/mL')

plt.grid(True)

plt.hold(False)



# super slow data:

datadir = '/home/kedmond/Dropbox/ZimmViscometer/Data/BTF2/'

files = glob.glob(datadir+'*stir*.txt')

data1 = np.genfromtxt(files[0], delimiter = ',', skiprows=2)
data2 = np.genfromtxt(files[1], delimiter = ',', skiprows=2)
data3 = np.genfromtxt(files[2], delimiter = ',', skiprows=2)

wm = np.array(data1[:, 1]) * 180/np.pi

wr1 = np.array(data1[:, 2]) * 180/np.pi
wr2 = np.array(data2[:, 2]) * 180/np.pi
wr3 = np.array(data3[:, 2]) * 180/np.pi

# 2 and 3 are one set that were split b/c of technical issues w/ the serial port:
wr4 = np.append(wr2, wr3)


fig = plt.figure(facecolor='white')

plt.loglog(wm, wr1, '-o')
plt.hold(True)
plt.loglog(wm, wr4, '-o')

#plt.loglog(wmcal, wr1cal, '-o')
#plt.loglog(wmcal, wr2cal, '-o')
#plt.loglog(wmcal, wr3cal, '-o')

ranges = np.array([.1, 1000])

plt.loglog(ranges, ranges*1.5, '--')

plt.xlabel(r'$\omega_m (deg./sec.)$', fontsize=20)

plt.ylabel(r'$\omega_r (deg./sec.)$', fontsize=20)

plt.title('0.8 g/L PEO, 1.051 g/mL')

plt.grid(True)

plt.hold(False)




fig = plt.figure(facecolor='white')

plt.loglog(wm-wr1, wr1, '-o')
plt.hold(True)
plt.loglog(wm-wr4, wr4, '-o')

#plt.loglog(wmcal-wr1cal, wr1cal, '-o')
#plt.loglog(wmcal-wr2cal, wr2cal, '-o')
#plt.loglog(wmcal-wr3cal, wr3cal, '-o')


ranges = np.array([.1, 1000])

plt.loglog(ranges, ranges*1.5, '--')

plt.xlabel(r'$\omega_m-\omega_r (deg./sec.)$', fontsize=20)

plt.ylabel(r'$\omega_r (deg./sec.)$', fontsize=20)

plt.title('0.8 g/L PEO, 1.051 g/mL')

plt.grid(True)

plt.hold(False)


files = glob.glob(datadir+'*.txt')

files.sort()

data1 = np.genfromtxt(files[0], delimiter = ',', skiprows=2)
data2 = np.genfromtxt(files[1], delimiter = ',', skiprows=2)
data3 = np.genfromtxt(files[2], delimiter = ',', skiprows=2)
data4 = np.genfromtxt(files[3], delimiter = ',', skiprows=2)
data5 = np.genfromtxt(files[4], delimiter = ',', skiprows=2)
data6 = np.genfromtxt(files[5], delimiter = ',', skiprows=2)

wr1 = np.array(data1[:, 2]) * 180/np.pi
wr2 = np.array(data2[:, 2]) * 180/np.pi
wr3 = np.array(data3[:, 2]) * 180/np.pi
wr4 = np.array(data4[:, 2]) * 180/np.pi
wr5 = np.array(data5[:, 2]) * 180/np.pi
wr6 = np.array(data6[:, 2]) * 180/np.pi

wm1 = np.array(data1[:, 1]) * 180/np.pi
wm2 = np.array(data2[:, 1]) * 180/np.pi
wm3 = np.array(data3[:, 1]) * 180/np.pi
wm4 = np.array(data4[:, 1]) * 180/np.pi
wm5 = np.array(data5[:, 1]) * 180/np.pi
wm6 = np.array(data6[:, 1]) * 180/np.pi

wr = (wr1 + wr2 + wr3[0:6])/3.0
wr = np.append(wr, wr3[6])
wrl = (wr4 + np.append(wr5, wr6))/2




dirname = '/home/kedmond/Dropbox/ZimmViscometer/Data/Calibration/'

files = glob.glob(dirname+'*0.8*.txt')

calib1 = np.genfromtxt(files[0], delimiter = ',', skiprows=2)
calib2 = np.genfromtxt(files[1], delimiter = ',', skiprows=2)
calib3 = np.genfromtxt(files[2], delimiter = ',', skiprows=2)

wm = np.array(calib1[:, 1]) * 180/np.pi

wr1 = np.array(calib1[:, 2]) * 180/np.pi
wr2 = np.array(calib2[:, 2]) * 180/np.pi
wr3 = np.array(calib3[:, 2]) * 180/np.pi



fig = plt.figure(facecolor='white')

ranges = np.array([.1, 1000])

plt.loglog(ranges, ranges, '--')

plt.hold(True)
plt.loglog(wm3, wr, '-o')
plt.loglog(wm4, wrl, '-o')

plt.loglog(ranges, ranges/2., '--')


plt.xlabel(r'$\omega_m (deg./sec.)$', fontsize=20)

plt.ylabel(r'$\omega_r (deg./sec.)$', fontsize=20)

plt.title('0.8 g/L PEO, 1.051 g/mL')

plt.grid(True)

plt.hold(False)

