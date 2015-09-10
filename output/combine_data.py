import numpy as np
import matplotlib.pyplot as plt

nobs = 13
nobspit = 1
date,mag,merr = [],[],[]
for iob in range(nobs):
    data = np.loadtxt('ob150448-%d.dat'%(iob+1))
    if iob == nobspit:
        np.savetxt('ob150448-space.dat',np.vstack([data[:,0],data[:,3],data[:,4]]).T,fmt='%f')
        continue
    date.extend(data[:,0])
    mag.extend(data[:,3])
    merr.extend(data[:,4])

date = np.array(date)
mag  = np.array(mag)
merr = np.array(merr)
order = np.argsort(date)
date = date[order]
mag = mag[order]
merr = merr[order]

ax = plt.subplot(111)
plt.errorbar(date,mag,yerr=merr,linestyle='-',color='k',marker='o')
np.savetxt('ob150448-ground.dat',np.vstack([date,mag,merr]).T,fmt='%f')
ax.invert_yaxis()
plt.show()
