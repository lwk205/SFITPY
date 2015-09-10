import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator,FormatStrFormatter
from PyAstronomy import pyasl
import matplotlib.cm as cm

#eventname = raw_input('Please give the eventname: ')
#eventname = eventname.upper()
#opts = open('figctrls.dat','r')
#for line in opts:
#    temp = line.split()
#    if temp[0] == eventname:
#        teams = temp[1].split('-')
#        cstrs = temp[2].split('-')
#        limits = np.array(temp[3:]).astype(float)

eventname = 'ob150845'
#teams = ['OGLE','Spitzer','CTIO I']
teams = ['OGLE','Spitzer','CTIO I','MOA','Danish Z','Danish V']
cstrs = cm.rainbow(np.linspace(0,1,len(teams)))
limits = [7150,7260,18.3,14.8,0.2]

fig = plt.figure(figsize=(8.5,11))
CJD = pyasl.get_juldate()-2450000 ## Current Julian Date
ax1 = plt.subplot(211)
ax2 = plt.subplot(212,sharex=ax1)
ax1.set_title(eventname)
#ax1.axvline(CJD,linestyle='--',color='k')
#ax2.axvline(CJD,linestyle='--',color='k')

eventname = eventname.lower()
for i in range(len(teams)):
    data = np.loadtxt('output/%s-%d.dat'%(eventname,i+1))
    ax1.errorbar(data[:,0],data[:,3],yerr=data[:,4],marker='o',color=cstrs[i],markersize=5,markeredgecolor=cstrs[i],markerfacecolor='none',linestyle='none',label=teams[i])
    ax2.errorbar(data[:,0],data[:,3]-data[:,5],yerr=data[:,4],marker='o',color=cstrs[i],markersize=5,markeredgecolor=cstrs[i],markerfacecolor='none',linestyle='none')
    if teams[i] in ['OGLE','Spitzer','CTIO H']:
        model = np.loadtxt('output/%s-%d-fin.dat'%(eventname,i+1))
        ax1.plot(model[:,0],model[:,3],color=cstrs[i])
        ax2.axhline(0,color='k')

## bin/plant model ##
#mod = np.loadtxt('ob150966-Bozza.mod')
#ax1.plot(mod[:,0],mod[:,1],color='k')
#ax2.plot(mod[:,0],mod[:,1],color='k')

################
tmin,tmax,Imax,Imin,merr = limits

ax1.text(CJD,(Imax+Imin)/2.,'current JD',rotation=90,va='center')
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#ax1.set_xlim(tmin,tmax)
#ax2.set_xlim(tmin,tmax)
#ax2.set_ylim(-0.2,0.2)
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_minor_locator(AutoMinorLocator(5))

ax1.legend(loc=0,numpoints=1)
ax1.axis([tmin,tmax,Imax,Imin])
ax2.axis([tmin,tmax,-merr,merr])
#ax1.invert_yaxis()
#ax2.invert_yaxis()
ax1.set_xlabel('HJD')
ax1.set_ylabel('I (OGLE)')
ax2.set_xlabel('HJD')
ax2.set_ylabel('Residuals')

plt.savefig('output/%s-lc.png'%eventname)
os.system('convert output/%s-lc.png output/%s-lc.pdf'%(eventname,eventname))
plt.show()
