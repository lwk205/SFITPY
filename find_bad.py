import numpy as np
import os

nobs = 6
eventname = 'ob150845'
chi2_crit = 20 #eventually, this threshold should depend on the number of data points

for iob in range(nobs):
    data = np.loadtxt('output/%s-%d.dat'%(eventname,iob+1))
    chi2 = (data[:,3]-data[:,5])**2/data[:,4]**2
    badpts = (chi2>chi2_crit)+(data[:,4]>0.5)
    bad = data[badpts,:]
    bad = np.vstack([bad.T,chi2[badpts]]).T
    if len(bad) != 0:
        np.savetxt('baddata/%s-%d.dat'%(eventname,iob+1),bad,fmt='%f')
