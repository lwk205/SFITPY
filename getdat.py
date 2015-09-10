#!/usr/bin/env python
import numpy as np

def getdat_sub(datfile,errfac,SPIT,tbegin,MJD):
### Here I assume the bad data have been removed ###
### TBA: identify and remove bad data ###
### TBA: MOA data need treated differently ###
    print 'Read in data file: ',datfile,SPIT,MJD
    data = np.loadtxt(datfile)
    date = data[:,0]
    mag  = data[:,1]
    merr = data[:,2]*errfac
    if date[0] > 2450000:
        date -= 2450000.
    good = date>tbegin
    date = date[good]
    mag  = mag[good]
    merr = merr[good]
#    anomaly = (date>7189.82)*(date<7190.91)
#    merr[anomaly] /= 4.
    if SPIT == True and MJD == True:
        print 'Spitzer use MJD'
        date += 0.5
    flux = 10.**(0.4*(18.-mag))
    ferr = merr*flux*np.log(10.)/2.5
    return np.vstack([date,flux,ferr,mag,merr])

def getalldat(datfiles,errfacs,isspit,parmref):
    tbegin,talert,t0par,t0base,MJD = parmref
    data = []
    print '---------Begin reading data...---------'
    for i in range(len(datfiles)):
        idata = getdat_sub(datfiles[i],errfacs[i],isspit[i],tbegin,MJD)
        data.append(idata)
    print '---------Reading data finished----------'
    return data

def rmbaddat(eventname,datfiles,data):
    nobs = len(datfiles)
    data_new = []
    for iob in range(nobs):
        try:
            baddata = np.loadtxt('baddata/%s-%d.dat'%(eventname,iob+1))
        except IOError:
            print 'no bad data has been found: nob=%d'%(iob+1)
            data_new.append(data[iob])
            continue
        date = data[iob][0]
        loc = []
        shape = baddata.shape
        if len(shape) >1:
            baddate = baddata[:,0]
        else:
            baddate = [baddata[0]]
        for ibad in baddate:
            iloc = np.abs(date-ibad)<1e-6
            if len(np.ones_like(iloc)[iloc]) != 1:
                print 'wrong baddata file'
                break
            loc.append(iloc)
        loc = np.sum(np.array(loc),axis=0).astype(bool)
        good = np.ones_like(loc).astype(bool)*(~loc)
        data_new.append((data[iob].T)[good].T)
    print '---------Bad data removed-----------'
    return np.array(data_new)

def convert_to_ogle_dat(flux,ferr,parmflx,parmflx_ogle,use_mag=True):
    fs,fb = parmflx
    fsogle,fbogle = parmflx_ogle
    flux_ogle = (flux-fb)/fs*fsogle+fbogle
    ferr_ogle = ferr/fs*fsogle
    if use_mag == False:
        return flux_ogle,ferr_ogle
    mag_ogle  = 18.-2.5*np.log10(flux_ogle)
    merr_ogle = ferr_ogle*2.5/np.log(10)/flux_ogle
    return mag_ogle,merr_ogle

def convert_to_ogle_mod(flux,parmflx,parmflx_ogle,use_mag=True):
    fs,fb = parmflx
    fsogle,fbogle = parmflx_ogle
    flux_ogle = (flux-fb)/fs*fsogle+fbogle
    if use_mag == False:
        return flux_ogle
    mag_ogle  = 18.-2.5*np.log10(flux_ogle)
    return mag_ogle

#def convert_flx_to_ogle_flx(flux,parmflx,parmflx_ogle):
#    fs,fb = parmflx
#    fsogle,fbogle = parmflx_ogle
#    flux_ogle = (flux-fb)/fs*fsogle+fbogle
#    return flux_ogle
