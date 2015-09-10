#!/usr/bin/env python
import numpy as np
import os
import getpar
import getbmag
import getmag_wei
import caustic_double_precision as gencau
import time

### Spitzer positions on the sky ###
spitz_traj = np.loadtxt('andy/spitz.13-15')
hjds = spitz_traj[:,0] + 6830.
ras  = spitz_traj[:,1]
decs = spitz_traj[:,2]
diss = spitz_traj[:,3]
nspit= len(hjds)
### File used for the finite source effect ###
b0b1 = np.loadtxt('andy/b0b1.dat')
b0tab = b0b1[:,1]
b1tab = b0b1[:,2]
db0tab= b0b1[:,3]
db1tab= b0b1[:,4]
#######################

def getb0p(z):
    b0 = np.zeros_like(z)
    b1 = np.zeros_like(z)
    db0= np.zeros_like(z)
    db1= np.zeros_like(z)
    clo = z<0.001
    med = (z>=0.001)*(z<10)
    far = z>=10.
    ## close part ##
    b0[clo] = 2*z[clo]
    b1[clo] = -5/14.*z[clo]
    db0[clo]= 2.
    db1[clo]= -5/14.
    ## far part ##
    b0[far] = 1+0.125/z[far]**2
    db0[far]= -0.25/z[far]**3
    b1[far] = 0.025/z[far]**2
    db1[far]= -0.05/z[far]**3
    ## intermediate part; this is different from the Fortran code ##
    k = (1000*z[med]).astype(int) #begins from 1; int=floor in Python, int=round in Fortran
    w1 = k+1-1000*z[med]
    w2 = 1.-w1
    b0[med] = w1* b0tab[k-1] + w2* b0tab[k]
    db0[med]= w1*db0tab[k-1] + w2*db0tab[k]
    b1[med] = w1* b1tab[k-1] + w2* b1tab[k]
    db1[med]= w1*db1tab[k-1] + w2*db1tab[k]
    ####
    return b0,b1,db0,db1

def lensMotion(te,piex,piey):
    ## If the lens undergoes circular motion ##
    period = 5. #period of the motion, days
    mratio = 0.001 #secondary-to-primary mass ratio
    murel  = 5.  #relative proper motion, mas/yr
    thetae = te*murel/365.25 # mas
    pie    = np.sqrt(piex**2+piey**2)
    pirel  = thetae*pie # mas
    Ml     = thetae/pie/8.14 #M_sun
    Dl     = 1./(pirel+0.125) # kpc
    sma    = mratio*(period**2*thetae/pie)**(1/3.)*0.01 # AU
    amp    = sma/thetae/Dl
    print '---Lens information---'
    print 'M_L (M_sun) = %.2f'%Ml
    print 'Period (days) = %.2f '%period
    print 'Secondary: q,s = %e,%.3f'%(mratio,sma/mratio/(Dl*thetae))
    print 'Mu_rel (mas/yr) = %.2f'%murel
    print 'Theta_E (mas) = %.2f'%thetae
    print 'D_L (kpc) = %.2f'%Dl
    print 'sma (AU) = %e'%sma
    print 'amp (theta_E) = %e'%amp
    print '----------------------'
    return period,amp

def getqnqe(dates,alpha,delta,t0par,qlat,qlong,SPITZ):
    qns,qes = [],[]
    qnqe = []
    for t in dates:
        qn,qe = getpar.geta(t,alpha,delta,t0par)
        if SPITZ == True:
            qnp,qep,qr = getpar.gets(t,hjds,ras,decs,diss,alpha,delta)
        else:
            qnp,qep = getpar.gett(t,qlat,qlong,alpha,delta)
        qn += qnp
        qe += qep
        qnqe.append([qn,qe])
    qnqe = np.array(qnqe).T
    return qnqe

def gettraj(t,t0,u0,te,piex,piey,alpha,delta,t0par,qlat,qlong,hjds,ras,decs,diss,SPITZ):
    tau = (t-t0)/te
    beta = u0
    qn,qe = getpar.geta(t,alpha,delta,t0par)
    if SPITZ == True:
        qnp,qep,qr = getpar.gets(t,hjds,ras,decs,diss,alpha,delta)
    else:
        qnp,qep = getpar.gett(t,qlat,qlong,alpha,delta)
    qn += qnp
    qe += qep
    dtau = piex*qn + piey*qe
    dbeta= piex*qe - piey*qn
    return tau+dtau,beta+dbeta

def getslc(alpha,delta,qlat,qlong,dates,parmfit,parmflx,Gamma,t0par,SPITZ,LENS,use_fsfb,qne=[]):
    ## piex == pien, piey = piee ##
    t0,u0,te,rho,piex,piey = parmfit
    fs,fb = parmflx
    traj = []
#    start = time.clock()
    if qne == []:
        for i in range(len(dates)):
            t = dates[i]
            taup,betap = gettraj(t,t0,u0,te,piex,piey,alpha,delta,t0par,qlat,qlong,hjds,ras,decs,diss,SPITZ)
            traj.append([taup,betap])
        traj = np.array(traj)
    else:
        tau = (dates-t0)/te
        beta= u0
        dtau = piex*qne[0]+piey*qne[1]
        dbeta= piex*qne[1]-piey*qne[0]
        traj = np.vstack([tau+dtau,beta+dbeta]).T
#    step1 = time.clock()
    x = np.sqrt(traj[:,0]**2+traj[:,1]**2)
    z = x/rho
    plps = (x**2+2)/x/np.sqrt(x**2+4)
    b0,b1,db0,db1 = getb0p(z)
    ## The use of this numerical method to compute the FS effect won't work if rho>0.1 ##
    if use_fsfb == True:
        famp = fs*plps*(b0-Gamma*b1)+fb
    else:
        famp = plps*(b0-Gamma*b1)
#    step2 = time.clock()
#    print step1-start,step2-step1
    return traj,famp


def getbinlc(alpha,delta,qlat,qlong,dates,parmfit,parmflx,Gamma,t0par,SPITZ,use_fsfb,caustic=[]):
    t0,u0,te,rho,piex,piey,theta,q,s = parmfit
    fs,fb = parmflx
    rho2 = rho**2
    theta = theta/180.*np.pi  #convert degree to radian
    if caustic == []: #if no caustic found, calculate it first
        ## lens positions centered on the center of magnification ##
        if s<1:
            offset_magcen = q/(1+q)*s
        else:
            offset_magcen = q/(1+q)/s
        x1 = -offset_magcen
        x2 = s-offset_magcen
        zcau,zcr = gencau.getCaustic(np.array([1.,q]),np.array([x1,x2]),np.array([0.,0.]),1000)
        xcau = np.real(zcau)
        ycau = np.imag(zcau)
    else: #otherwise use the available caustic
        xcau = caustic[0]
        ycau = caustic[1]
##############################
    if s < 1.:
        offset_x = s*(-0.5+1./(1+q))
    else:
        offset_x = s/2.-q/(1+q)/s
#    offset_x = s/2.
    magmap = []
    for i in range(len(dates)):
        taup,betap = gettraj(dates[i],t0,u0,te,piex,piey,alpha,delta,t0par,qlat,qlong,hjds,ras,decs,diss,SPITZ)
        xcm = taup*np.cos(theta) + betap*np.sin(theta)
        ycm =-taup*np.sin(theta) + betap*np.cos(theta)
        xs = xcm-offset_x
        ys = ycm
        dist2 = min((xcau-xcm)**2+(ycau-ycm)**2)
#        dist2 = min((xcau-xs)**2+(ycau-ys)**2)
        if dist2 < 4*rho2:
        ### Jin code ###
            magbps,magnew,errorflag = getbmag.getmag_jin(xs,ys,s,q,rho,Gamma)
            magnew = -1
            if magnew < 1.:
                print 'looplinking failed, switch to mapmaking'
                magnew,errorflag = getmag_wei.getmag_wei(magbps,xs,ys,s,q,rho,Gamma)
        elif dist2 < 16*rho2:
        ### hex approx ###
            magnew = getbmag.taylor_2(xs,ys,q,s,Gamma,rho,4)
        elif dist2 < 100*rho2:
        ### quad approx ###
            magnew = getbmag.taylor_2(xs,ys,q,s,Gamma,rho,2)
        else:
        ### mon approx ###
            magnew = getbmag.taylor_2(xs,ys,q,s,Gamma,rho,0)
        if use_fsfb == True:
            magnew = magnew*fs+fb
        magmap.append([xcm,ycm,magnew])
    return np.array(magmap),xcau,ycau

