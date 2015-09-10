#!/usr/bin/env python
import numpy as np
import genlc
import emcee
import getdat
import getparm
from scipy.optimize import fmin
import time
#import matplotlib.pyplot as plt

def initialization(parms_names,isspitz):
    parmfit,namefit,parmflx,nameflx,parmref = parms_names
    ## switch free_fsfb to False if use getchi2 !!!
    free_fsfb = False
    if True in isspitz:
        namefix = np.array(['rho'])
    else: #only ground-based data: set parallax fixed
        namefix = np.array(['rho','pien','piee'])
#    namefix = np.array(['rho'])
    if free_fsfb == False:
        namefix = np.hstack([namefix,nameflx.flatten()])
    nameall = np.hstack([namefit,nameflx.flatten()])
    parmall = np.hstack([parmfit,parmflx.flatten()])
    kfree = np.ones_like(nameall).astype(bool)
    for iname in namefix:
        kfree[nameall==iname] = False
    parmfree = parmall[kfree]
    namefree = nameall[kfree]
    parmfix = parmall[~kfree]
    print 'Free parameters:'
    print namefree
    print parmfree
    print 'Fixed parameters:'
    print namefix
    print parmfix
    print '------ Initialization finished ------'
    return parmfree,namefree,parmfix,namefix

def getparmfree(parmall,nameall,kfix):
    parmfree = parmall[kfix==False]
    namefree = nameall[kfix==False]
    parmfix  = parmall[kfix==True]
    namefix  = nameall[kfix==True]
    return parmfree,parmfix,namefree,namefix

def getparmfit(parmfree,parmfix,namefree,namefix):
    parmall = np.hstack([parmfree,parmfix])
    nameall = np.hstack([namefree,namefix])
    namefit = np.array(['t0','u0','te','rho','pien','piee'])
    parmfit = [parmall[nameall==i] for i in namefit]
    parmflx = [parmall[nameall==i] for i in nameall if not (i in namefit)]
    parmfit = np.array(parmfit).flatten()
    parmflx = np.array(parmflx).reshape((-1,2))
#    print parmfit
#    print parmflx
    return parmfit,parmflx

### 'getfsfb' and 'getchi2' works ###
def getfsfb(iflux,iferr,iamp):
    a = np.ones(2)
    for iloop in range(2):
        chi2 = 0.
        d = np.zeros(2)
        b = np.zeros((2,2))
        f = np.zeros(2)
        for k in range(len(iamp)):
            y = iflux[k]
            sig2 = iferr[k]**2
            f[0] = 1.
            f[1] = iamp[k]
            ymod = 0.
            for i in range(2):
                ymod += a[i]*f[i]
                d[i] += y*f[i]/sig2
                for j in range(2):
                    b[i,j] += f[i]*f[j]/sig2
            chi2add = (y-ymod)**2/sig2
            chi2 += chi2add
        if iloop == 1:
            continue
        c = np.linalg.inv(b)
        for i in range(2):
            a[i] = 0.
            for j in range(2):
                a[i] += c[i,j]*d[j]
    fb = a[0]
    fs = a[1]
    fberr = np.sqrt(c[0,0])
    fserr = np.sqrt(c[1,1])
    return chi2,fs,fb,fserr,fberr

#def getchi2(parmfree,parmfix,namefree,namefix,data,nobs,eventinfo,t0par,isspit,ischi2,get_chi2_obs=False):
#    date,flux,ferr = data
#    alpha,delta,oblat,oblong,obgamma = eventinfo
#    parmfit,parmflx = getparmfit(parmfree,parmfix,namefree,namefix)
#    chi2 = []
#    fsfb = []
#    ## apply the color constraint (2-sigma limit) ##
#    I_minus_L = [-1.2,-0.7]
#    fsogle = parmflx[0,0]
#    fsspit = parmflx[1,0]
#    I_minus_L_trial = 2.5*np.log10(fsspit/fsogle)
#    if (I_minus_L_trial-I_minus_L[0])*(I_minus_L_trial-I_minus_L[1]) <= 0:
#        chi2_color = 0.
#    else:
#        chi2_color = 0.#np.inf
#    for iob in range(nobs):
#        itraj,iamp = genlc.getslc(alpha,delta,oblat[iob],oblong[iob],date[iob],parmfit,parmflx[iob],obgamma[iob],t0par,isspit[iob],False,False)
#        ichi2,fs,fb,fserr,fberr = getfsfb(flux[iob],ferr[iob],iamp)
#        ## negative blending (fb<-0.2) should be avoided ##
#        if parmflx[iob][1] < -0.2:
#            chi2_blend = 0.#np.inf
#        else:
#            chi2_blend = 0.
#        chi2.append(ichi2+chi2_blend)
#        fsfb.append([fs,fb,fserr,fberr])
#    if get_chi2_obs == True:
#        return np.array(chi2),np.array(fsfb)
#    chi2tot = np.sum(chi2)+chi2_color
#    if ischi2 == True:
#        return chi2tot
#    else:
#        return -0.5*chi2tot

def getchi2(parmfree,parmfix,namefree,namefix,data,eventinfo,t0par,isspitz,ischi2,use_color,I_minus_L,qnqe,get_chi2_obs=False):
    alpha,delta,oblats,oblongs,obgammas = eventinfo
    parmfit,parmflx = getparmfit(parmfree,parmfix,namefree,namefix)
    dof,chi2 = [],[]
    ######################
    for iob in range(len(data)):
        date,flux,ferr = data[iob][:3]
        qnqe_iob = qnqe[iob]
        itraj,ifmod = genlc.getslc(alpha,delta,oblats[iob],oblongs[iob],date,parmfit,parmflx[iob],obgammas[iob],t0par,isspitz[iob],False,False,qne=qnqe_iob)
        ichi2,fs,fb,fserr,fberr = getfsfb(flux,ferr,ifmod)
        chi2.append(ichi2)
        parmflx[iob][0] = fs
        parmflx[iob][1] = fb
        dof.append(len(date))
    ## apply the color constraint (2-sigma limit) ##
    chi2_color = 0.
    if use_color==True:
        fsogle = parmflx[0][0]
        fsspit = parmflx[isspitz==True][0][0]
        I_minus_L_trial = 2.5*np.log10(fsspit/fsogle)
        if (I_minus_L_trial-I_minus_L[0])*(I_minus_L_trial-I_minus_L[1]) > 0:
            chi2_color = np.inf
    ## negative blending (fb<-0.2) should be avoided: only apply to OGLE ##
    if parmflx[0][1] < -0.2:
        chi2_blend = np.inf
    else:
        chi2_blend = 0.
    if get_chi2_obs == True:
        return np.array(chi2),np.array(dof),parmflx
    chi2tot = np.sum(chi2)
    if ischi2 == True:
#        print chi2tot+chi2_color+chi2_blend
        return chi2tot+chi2_color+chi2_blend
    else:
        return -0.5*chi2tot

def lnprob_func(parmfree,parmfix,namefree,namefix,data,eventinfo,t0par,isspitz,ischi2,use_color,I_minus_L,qnqe,get_chi2_obs=False):
    alpha,delta,oblats,oblongs,obgammas = eventinfo
    parmfit,parmflx = getparmfit(parmfree,parmfix,namefree,namefix)
#    namefit = np.array(['t0','u0','te','rho','pien','piee'])
#    u0 = parmfit[namefit=='u0'][0]
    dof,chi2 = [],[]
    ## apply the color constraint (2-sigma limit) ##
    chi2_color = 0.
    if use_color==True:
        fsogle = parmflx[0][0]
        fsspit = parmflx[isspitz==True][0][0]
        I_minus_L_trial = 2.5*np.log10(fsspit/fsogle)
        if (I_minus_L_trial-I_minus_L[0])*(I_minus_L_trial-I_minus_L[1]) > 0:
            chi2_color = np.inf
    ## negative blending (fb<-0.2) should be avoided: only apply to OGLE ##
    if parmflx[0][1] < -0.2:
        chi2_blend = np.inf
    else:
        chi2_blend = 0.
    ######################
    for iob in range(len(data)):
        date,flux,ferr = data[iob][:3]
        qnqe_iob = qnqe[iob]
        itraj,ifmod = genlc.getslc(alpha,delta,oblats[iob],oblongs[iob],date,parmfit,parmflx[iob],obgammas[iob],t0par,isspitz[iob],False,True,qne=qnqe_iob)
        chi2.append(np.sum((flux-ifmod)**2/ferr**2))
        dof.append(len(date))
    if get_chi2_obs == True:
        return np.array(chi2),np.array(dof),parmflx
    chi2tot = np.sum(chi2)
    if ischi2 == True:
#        print chi2tot+chi2_color+chi2_blend
        return chi2tot+chi2_color+chi2_blend
    else:
        return -0.5*chi2tot

##########################################
##### Fitting program based on EMCEE #####
##########################################
def mcmc(prob_func,eventname,namefree,parmfree,namefix,parmfix,data,eventinfo,parmref,isspitz,use_color=False,I_minus_L=[]):
    start = time.clock()
    tbegin,talert,t0par,t0base,mjd = parmref
    alpha,delta,oblats,oblongs,obgammas = eventinfo
    ## compute (qn,qe) once and for all ##
    qnqe = []
    for iob in range(len(data)):
        date = data[iob][0]
        qnqe_iob = genlc.getqnqe(date,alpha,delta,t0par,oblats[iob],oblongs[iob],isspitz[iob])
        qnqe.append(qnqe_iob)
### Get chi2 based on initial parameters ###
#    chi2_obs,dof,fsfb = prob_func(parmfree,parmfix,namefree,namefix,data,eventinfo,t0par,isspitz,True,use_color,I_minus_L,qnqe,get_chi2_obs=True)
#    print chi2_obs,dof
    ndim = len(parmfree)
    nwalkers = 30
    pos = [parmfree+1e-3*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers,ndim,prob_func,args=(parmfix,namefree,namefix,data,eventinfo,t0par,isspitz,False,use_color,I_minus_L,qnqe),threads=1)
    pos,lnprob,rstate = sampler.run_mcmc(pos,100)  ## burn in
    sampler.reset()
    sampler.run_mcmc(pos,200)  ## MCMC sampling
    ## save the MCMC chain and other useful information ##
    chain = sampler.chain.reshape((-1,ndim),order='F')
    chi2s = -2*sampler.lnprobability.reshape(-1,order='F')
    headerstr = '%s '%namefree[0]
    for i in range(1,len(namefree)):
        headerstr += '%s '%namefree[i]
    np.savetxt('output/mcmc/chain-%s.dat'%eventname,chain,fmt='%f',header=headerstr)
    np.savetxt('output/mcmc/acceptance-%s'%eventname,sampler.acceptance_fraction,fmt='%f')
    np.savetxt('output/mcmc/chisq-%s.dat'%eventname,chi2s)
    ## find the best-fit and error bars ##
    parmbest = chain[chi2s==min(chi2s)][0]
    parmmean = np.median(chain,axis=0)
    parmerrs = np.std(chain,axis=0)
    print 'best parameters: '
    print namefree
    print ['%.4f'%iparm for iparm in parmbest]
    print ['%.4f'%iparm for iparm in parmmean]
    print ['%.4f'%iparm for iparm in parmerrs]
    chi2_obs,dof,fsfb = prob_func(parmbest,parmfix,namefree,namefix,data,eventinfo,t0par,isspitz,True,use_color,I_minus_L,qnqe,get_chi2_obs=True)
    print 'MCMC finished (s): ',time.clock()-start,' s'
    chi2min = np.sum(chi2_obs)
#    return parmbest,parmfix,np.hstack([chi2min,chi2_obs]),dof
    parmfit,parmflx = getparmfit(parmbest,parmfix,namefree,namefix)
    parmflx = fsfb
    return parmfit,parmflx,np.hstack([chi2min,chi2_obs]),dof

#############################################
##### Fitting program based on downhill #####
#############################################
def downhill(prob_func,namefree,parmfree,namefix,parmfix,data,eventinfo,parmref,isspitz,use_color=False,I_minus_L=[]):
    start = time.clock()
    tbegin,talert,t0par,t0base,mjd = parmref
    alpha,delta,oblats,oblongs,obgammas = eventinfo
    ## compute (qn,qe) once and for all ##
    qnqe = []
    for iob in range(len(data)):
        date = data[iob][0]
        qnqe_iob = genlc.getqnqe(date,alpha,delta,t0par,oblats[iob],oblongs[iob],isspitz[iob])
        qnqe.append(qnqe_iob)
### Get chi2 based on initial parameters ###
#    chi2_obs,dof,fsfb = prob_func(parmfree,parmfix,namefree,namefix,data,eventinfo,t0par,isspitz,True,use_color,I_minus_L,qnqe,get_chi2_obs=True)
#    print chi2_obs,dof
#    parmbest = fmin(prob_func,parmfree,args=(parmfix,namefree,namefix,data,eventinfo,t0par,isspitz,True,use_color,I_minus_L))
    parmbest,chi2min,iter,funcalls,warnflag,allevcs = fmin(prob_func,parmfree,args=(parmfix,namefree,namefix,data,eventinfo,t0par,isspitz,True,use_color,I_minus_L,qnqe),full_output=True,retall=True,maxiter=3000,maxfun=10000)
    print 'best parameters: '
    print namefree
    print ['%.3f'%iparm for iparm in parmbest]
    chi2_obs,dof,fsfb = prob_func(parmbest,parmfix,namefree,namefix,data,eventinfo,t0par,isspitz,True,use_color,I_minus_L,qnqe,get_chi2_obs=True)
    print 'Downhill minimization finished (s): ',time.clock()-start,' s, ',iter,' steps'
    chi2min = np.sum(chi2_obs)
    parmfit,parmflx = getparmfit(parmbest,parmfix,namefree,namefix)
    parmflx = fsfb
    return parmfit,parmflx,np.hstack([chi2min,chi2_obs]),dof
#    return parmbest,parmfix,np.hstack([chi2min,chi2_obs]),dof

def massfit_obj(q,s,thetas,lcs,lcs_full,eventinfo,parms_names,isspitz,use_color,I_minus_L,outputlc):
    alpha,delta,oblats,oblongs,obgammas = eventinfo
    parmfit,namefit,parmflx,nameflx,parmref = parms_names
    tbegin,talert,t0par,t0base,mjd = parmref
    parmfree,namefree,parmfix,namefix = initialization(parms_names,isspitz)
    sinlcs,chi2all = [],[]
    nlc = len(lcs)
    bestfits = []
    for i in range(nlc):
        data = lcs[i]
        iparmfit,iparmflx,chi2,dof = downhill(lnprob_func,namefree,parmfree,namefix,parmfix,data,eventinfo,parmref,isspitz,use_color=use_color,I_minus_L=I_minus_L)
        ### find best-fit model ###
        if outputlc == True:
#            iparmfit,iparmflx = getparmfit(iparmfree,parmfix,namefree,namefix)
            for iob in range(len(data)):
                date,flux,ferr = data[iob][:3]
                itraj,ifmod = genlc.getslc(alpha,delta,oblats[iob],oblongs[iob],date,iparmfit,iparmflx[iob],obgammas[iob],t0par,isspitz[iob],False,True)
#                plt.errorbar(date,flux,yerr=ferr,linestyle='none')
#                plt.plot(date,ifmod)
#                plt.show()
                np.savetxt('output/models/fit-ilc=%d-iob=%d.dat'%((i+1),(iob+1)),np.vstack([lcs_full[i][iob],ifmod]),fmt='%f',header='fs=%f, fb=%f'%(iparmflx[iob][0],iparmflx[iob][1]))
        bestfits.append(np.hstack([q,s,thetas[i],chi2.reshape(-1),iparmfit,iparmflx]))
        print '##########',q,s,thetas[i],chi2[0],'############'
    return bestfits

def massfit_sub(q,s,thetas,lcs,lcs_full,eventinfo,parms_names,isspitz,use_color,I_minus_L,outputlc):
    alpha,delta,oblats,oblongs,obgammas = eventinfo
    parmfit,namefit,parmflx,nameflx,parmref = parms_names
    tbegin,talert,t0par,t0base,mjd = parmref
    ## find data released before alert date ##
    nlc = len(lcs)
    lcs_sub = []
    oblats_sub,oblongs_sub,obgammas_sub,isspitz_sub = [],[],[],[]
    parmflx_sub,nameflx_sub = [],[]
    for i in range(nlc):
        data = lcs[i]
        data_sub = []
        flag = False
        for iob in range(len(data)):
            date,flux,ferr = data[iob][:3]
            good = date<talert
            if len(date[good]) == 0:
                flag = True
                continue
            if i == 0:
                oblats_sub.append(oblats[iob])
                oblongs_sub.append(oblongs[iob])
                obgammas_sub.append(obgammas[iob])
                isspitz_sub.append(isspitz[iob])
                parmflx_sub.append(parmflx[iob])
                nameflx_sub.append(nameflx[iob])
            data_sub.append([date[good],flux[good],ferr[good]])
#        if flag == True:
#            continue
        lcs_sub.append(data_sub)
    oblats_sub = np.array(oblats_sub)
    oblongs_sub = np.array(oblongs_sub)
    obgammas_sub = np.array(obgammas_sub)
    isspitz_sub = np.array(isspitz_sub)
    parmflx_sub = np.array(parmflx_sub)
    nameflx_sub = np.array(nameflx_sub)
    eventinfo_sub = [alpha,delta,oblats_sub,oblongs_sub,obgammas_sub]
    parms_names_sub = [parmfit,namefit,parmflx_sub,nameflx_sub,parmref]
    if not (True in isspitz_sub):
        use_color_sub = False
    parmfree,namefree,parmfix,namefix = initialization(parms_names,isspitz)
    parmfree_sub,namefree_sub,parmfix_sub,namefix_sub = initialization(parms_names_sub,isspitz_sub)
    sinlcs,chi2all = [],[]
    bestfits = []
    nobs = len(lcs[0])
    nobs_sub = len(lcs_sub[0])
    for i in range(nlc):
        data_sub = lcs_sub[i]
        iparmfit_sub,iparmflx_sub,chi2_sub,dof_sub = downhill(lnprob_func,namefree_sub,parmfree_sub,namefix_sub,parmfix_sub,data_sub,eventinfo_sub,parmref,isspitz_sub,use_color=use_color_sub,I_minus_L=I_minus_L)
        if chi2_sub[0] > 10:
            chi2 = np.ones(1+nobs)
        else:
            data = lcs[i]
            iparmfit,iparmflx,chi2,dof = downhill(lnprob_func,namefree,parmfree,namefix,parmfix,data,eventinfo,parmref,isspitz,use_color=use_color,I_minus_L=I_minus_L)
#        ### find best-fit model ###
#        if outputlc == True:
#            iparmfit,iparmflx = getparmfit(iparmfree,parmfix,namefree,namefix)
#            for iob in range(len(data)):
#                date,flux,ferr = data[iob][:3]
#                itraj,ifmod = genlc.getslc(alpha,delta,oblats[iob],oblongs[iob],date,iparmfit,iparmflx[iob],obgammas[iob],t0par,isspitz[iob],False,True)
##                plt.errorbar(date,flux,yerr=ferr,linestyle='none')
##                plt.plot(date,ifmod)
##                plt.show()
#                np.savetxt('output/models/fit-ilc=%d-iob=%d.dat'%((i+1),(iob+1)),np.vstack([lcs_full[i][iob],ifmod]),fmt='%f',header='fs=%f, fb=%f'%(iparmflx[iob][0],iparmflx[iob][1]))
        bestfits.append(np.hstack([q,s,thetas[i],chi2.reshape(-1),chi2_sub.reshape(-1),iparmfit_sub,iparmflx_sub]))
#        bestfits.append(np.hstack([q,s,thetas[i],chi2.reshape(-1),chi2_sub.reshape(-1),iparmfree,parmfix,iparmfree_sub,parmfix_sub]))
        print '##########',q,s,thetas[i],chi2[0],'############'
    return bestfits
