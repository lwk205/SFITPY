import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import genlc
from matplotlib.ticker import MultipleLocator,AutoMinorLocator

def main():
    ## Event Information ##
    eventname = 'ob151346'
    alpha = (17. + 59/60. + 13.67/3600.)*15
    delta =-(28. + 56/60. + 11.8/3600.)
    fs,fb = 0.07
    Gamma = 0.60
    ## OGLE site ##
    qlong = -70.702
    qlat  = -29.0083
    ## planet/binary parameters ##
    t0 = 7190.74
    u0 = 0.18
    te = 25.12
    rho = 0.001
    theta = 80. # in degrees
    q = 0.0014
    s = 1.116
    Ibase = 19.4027
    fb_over_fs = -0.4134
#    ## Binary Parameters for ob150914 ##
#    t0 = 7172.14
#    u0 = 0.0243
#    te = 36.8418
#    rho = 0.0047
#    theta = 3.81532/np.pi*180
#    theta = 180-theta ## Bozza system
#    q = 0.369
#    s = 1.4718
#    ## Bozza gave Ibase and Fb/Fs ##
#    Ibase = 14.6422
#    fb_over_fs = 18.445
#    ## binary event ob150060 ##
#    t0 = 7148.109
#    u0 = 0.019
#    te = 71.8
#    rho = 0.00479
#    theta = 3.593/np.pi*180
#    q = 1.10
#    s = 3.40
#    Ibase = 16.6885
#    fb_over_fs = 0.3337
#    ## Planetary Event OB150051 ##
#    t0 = 7083.100
#    u0 = 0.221
#    te = 10.95
#    rho = 0.045
#    theta = 5.351*180/np.pi #Han
#    q = 5.94e-3
#    s = 0.939
#    Ibase = 16.629
#    fb_over_fs = -0.0216
    ###########
    fb_plus_fs = 10**(0.4*(18-Ibase))
    fs = fb_plus_fs/(1+fb_over_fs)
    fb = fb_over_fs*fs
    print fs,fb
    ###########
    fig = plt.figure(figsize=(8.5,11))
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    linestr = ['-','--','-.']
#    colors = ['k','#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','0.5']
#    pies = [-0.1,0,0.1]
    pieNs = [0.0193]
    pieEs = [-0.2949]
    colors = cm.rainbow(np.linspace(0,1,len(pieEs)*len(pieNs)))
    nmod = 0
#    ogle = np.loadtxt('data/ob150966.dat')
#    dates = ogle[:,0]-2450000
    dates = np.linspace(7100,7250,15000)
    for pieE in pieEs:
        for pieN in pieNs:
            piex = pieN
            piey = pieE
            parmfit = [t0,u0,te,rho,piex,piey,theta,q,s]
            parmflx = [fs,fb]
            t0par = t0
            lcg,xcau,ycau = genlc.getbinlc(alpha,delta,qlat,qlong,dates,parmfit,parmflx,Gamma,t0par,False,True)
            lcs,xcau,ycau = genlc.getbinlc(alpha,delta,qlat,qlong,dates,parmfit,parmflx,Gamma,t0par,True,True)
            magg = 18-2.5*np.log10(lcg[:,2])
            mags = 18-2.5*np.log10(lcs[:,2])
            if nmod == 0:
                ax1.plot(dates,magg,color='b',linestyle='--',label='OGLE',lw=2)
#                np.savetxt('ob150966-Bozza.mod',np.vstack([dates,magg]).T,fmt='%f')
                ax2.plot(dates,magg,color='b',linestyle='--',label='OGLE',lw=2)
                ax3.plot(xcau,ycau,marker='.',markersize=2,linestyle='none',color='k')
                ax3.plot(lcg[:,0],lcg[:,1],color='b')
            ax1.plot(dates,mags,color=colors[nmod],label='Spitzer,(%.1f,%.1f)'%(pieE,pieN))
            ax2.plot(dates,mags,color=colors[nmod],label='Spitzer,(%.1f,%.1f)'%(pieE,pieN))
            ax3.plot(lcs[:,0],lcs[:,1],color=colors[nmod])
            nmod += 1
    ax3.axis('equal')
    plt.suptitle(eventname,fontsize=14,fontweight='bold')
    ax1.legend(loc=0,fontsize=10)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax2.axis([7143,7156,15.3,11.7])
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax3.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax3.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.show()

main()
