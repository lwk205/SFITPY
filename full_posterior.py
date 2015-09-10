import numpy as np
import matplotlib.pyplot as plt
import triangle

#eventname = raw_input('Please give the event name (eg, ob140124): ')
eventname = 'ob150763'

chain = np.loadtxt('output/mcmc/chain-%s.dat'%eventname)
chain[:,0] -= int(chain[0,0])
chisq = np.loadtxt('output/mcmc/chisq-%s.dat'%eventname)

f = open('output/mcmc/chain-%s.dat'%eventname,'r')
strs = f.readline()
labelstrs = strs.split()[1:]

parms = [[min(chisq),0]]
for i in range(len(labelstrs)):
    parms.append([np.median(chain[:,i]),np.std(chain[:,i])])
np.savetxt('parms-%s.dat'%eventname,parms,fmt='%f',header=strs)

#labelstrs = [r'$t_{\rm ce}$ (HJD$^\prime-6851$)',
#r'$u_0$',
#r'$t_{\rm E}$ (days)',
#r'$\rho$ (10$^{-4}$)',
#r'$\pi_{{\rm E},N}$',
#r'$\pi_{{\rm E},E}$',
#r'$\alpha$ (deg)',
#r'$s$',
#r'$q$']
#triangle.corner(params.T,labels=labelstrs,plot_datapoints=False)

fig = triangle.corner(chain,labels=labelstrs,plot_datapoints=False)
fig.suptitle(eventname.upper())
plt.savefig('output/%s-posterior.pdf'%eventname)
plt.show()
