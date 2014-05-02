from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


'''
Created by Zach Forbes 23 Feb 2013

This program will calculate the viscosity of a fluid based on Zimm-Crothers viscometer
data. The first hard-coded example will be a background sample of NaBr in water with 
a solution of 1.75 um TPM hard spheres.

History:
    Updated by KVE, 3.06.2013

'''

def calcVisc(wr_0, wr, C, eta_0 = 0.0014):
    # viscosity of NaBr in H2O is 1.4 mPa*s at 20C (tabulated)
    # viscosity of H2O is 1 mPa s at 20C, 0.89 mPa s at 25C, 0.797 mPa s at 30C
    wr_0mean = np.mean(wr_0)
    wrmean = np.mean(wr)
    w = wr_0mean/wrmean
    
    viscosity = w*eta_0 + C*(w-1)
    
    return viscosity

def getData(directory, filename):
    data = np.genfromtxt(directory+filename, skip_header=2, delimiter=',', 
                         usecols = (0, 1, 2), dtype=('float64, float64, float64'))
    data.dtype.names = ('T', 'wm', 'wr')    # label the columns
    # sort in terms of temperature:
    data = np.sort(data, order='T')
    return data

def calculations(samplePath, sampleFile, conc, colors,
         backgroundDataPath='/media/data/Data/ZimmViscometer/KaZimmData/NaBr_H2O_Calibration/',
         backgroundDataFile='12-10-19-15-58-16_NaBr+H2O.txt'):

    background = getData(backgroundDataPath, backgroundDataFile)
    sample = getData(samplePath, sampleFile)
    backC, backTval, backTemps, backTnum = cValues(background, 'back')
    sampleC, sampleTval, sampleTemps, sampleTnum = cValues(sample, 'sample')

    visc = np.zeros(sampleTnum)
 
    for i in range(sampleTnum):
        wrSample = sample['wr'][i*sampleTval:(i+1)*sampleTval]
        wrUsable = np.zeros(12)
        for j in range(12):
            wrUsable[j] = wrSample[2+j]
        
        visc[i] = calcVisc(background['wr'], wrUsable, sampleC[i])
      
    return visc
 

def cValues(data,type, fit=True):    
    tval = len( np.where(data['T'] == data['T'][0])[0] )
    temps = np.unique(data['T'])
    tnum = len(temps)

    if fit:
        clist = np.zeros(tnum, dtype='float32')
    for i in range(0, tnum, 1):
        x = data['wm'][i*tval:(i+1)*tval] - data['wr'][i*tval:(i+1)*tval]
        y = data['wr'][i*tval:(i+1)*tval]

        if fit:
            if type == 'sample':
                clist[i] = findInflection(x,y)
            elif type == 'back':
                clist[i], _, _, _ = np.linalg.lstsq(x[:, np.newaxis], y)
            else: pass

    return clist, tval, temps, tnum


def findInflection(x, y, threshold=0.1):
    m, _, _, _ = np.linalg.lstsq(x[:, np.newaxis], y)
    yp = y - m * (x - np.mean(x))

    if np.std(yp) > 0.2:
        xp = x[yp == np.max(yp)]
        mf, _, _, _ = np.linalg.lstsq(x[x<xp][:, np.newaxis], y[x<xp])
    
        return mf
    else:
        return m

#viscosity_01 = main('/Users/Zach/Dropbox/Zimm viscometer/KaZimmData/TPM0.1-phi/',
#    '13-02-12-16-57-34_1.75um-TPM-0.1conc.txt','0.1','r')

files = ['13-02-21-12-25-15_1.75um-TPM-0.05conc-2.txt',
         '13-02-12-16-57-34_1.75um-TPM-0.1conc.txt',
         '13-01-29-11-03-31_1.75um-TPM.txt']
directories = ['/media/data/Data/ZimmViscometer/KaZimmData/TPM0.05-phi/',
               '/media/data/Data/ZimmViscometer/KaZimmData/TPM0.1-phi/',
               '/media/data/Data/ZimmViscometer/KaZimmData/KVE01-31B-TPM-0.2/']
concentrations = [0.05, 
                  0.1, 
                  0.2]
colors = ['r','g','b']


def main():
    plt.clf()
    
    v20 = np.zeros(len(files))
    v25 = np.zeros(len(files))
    v30 = np.zeros(len(files))
    
    for i in range(len(files)):
    
        v20[i] = calculations(directories[i],files[i], concentrations[i], colors[i])[0]
        v25[i] = calculations(directories[i],files[i], concentrations[i], colors[i])[1]
        v30[i] = calculations(directories[i],files[i], concentrations[i], colors[i])[2]
    
    v20 = v20.tolist()
    v25 = v25.tolist()
    v30 = v30.tolist()
    
    plt.plot(concentrations, v20, 'o', color=colors[0], label = '20$^\circ$ C')
    plt.plot(concentrations, v25, 'D', color=colors[1], label = '25$^\circ$ C')
    plt.plot(concentrations, v30, 's', color=colors[2], label = '30$^\circ$ C')

    plt.xlabel('$\phi$')
    plt.ylabel('$\eta_r$')
    plt.grid(True)
    plt.xlim(0, .25, .05)
    plt.ylim(0, 0.4)
    plt.legend(loc='lower left')
    plt.show()
main()