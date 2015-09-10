import numpy as np
import ConfigParser

def getlatlong(id_obs,datfile):
    if id_obs == 'ogle':
        qlong = -70.702
        qlat  = -29.0083
    elif id_obs == 'moa':
        qlong = +(170 + 27.9/60.)
        qlat  = -(43 + 59.2/60.)
    elif id_obs =='wise':
        qlong = +(34 + 45/60. + 44/3600.)
        qlat  = +(30 + 35/60. + 50/3600.)
    elif id_obs =='ct13':
        qlong = -70.815
        qlat  = -30.165
    elif id_obs == 'spitzer':
        qlong = 200
        qlat  = 100
    else:
        print 'lat-long not found: use OGLE for default'
        qlong = -70.702
        qlat = -29.0083
    return qlat,qlong

def getEventInfo(eventname):
    filename = 'event-%s.cfg'%eventname
    p = ConfigParser.RawConfigParser()
    p.read(filename)
    ename = p.get('Event Info','name')
    ra_str = p.get('Event Info','ra_j2000')
    dec_str = p.get('Event Info','dec_j2000')
    gamma_i = float(p.get('Source Info','gamma_i'))
    gamma_v = float(p.get('Source Info','gamma_v'))
    gamma_h = float(p.get('Source Info','gamma_h'))
    gamma_l = float(p.get('Source Info','gamma_l'))
    gamma_all = np.array([gamma_i,gamma_v,gamma_h,gamma_l])
    band_names= np.array(['I','V','H','L'])
    obsnames = p.get('Data Files','observatories').split()
    bands = p.get('Data Files','bands').split()
    nobs = len(obsnames)
    datfiles = []
    oblats,oblongs,obgammas = [],[],[]
    SPITZ = []
    for iob in range(nobs):
        datfile = p.get('Data Files',obsnames[iob])
        qlat,qlong = getlatlong(obsnames[iob],datfile)
        ispitz = False
        if obsnames[iob] == 'spitzer':
            ispitz = True
        gamma = gamma_all[band_names==bands[iob]][0]
        datfiles.append(datfile)
        oblats.append(qlat)
        oblongs.append(qlong)
        obgammas.append(gamma)
        SPITZ.append(ispitz)
    SPITZ = np.array(SPITZ)
    ra = np.array(ra_str.split(':')).astype(float)
    alpha = (ra[0]+ra[1]/60.+ra[2]/3600.)*15.
    dec = np.abs(np.array(dec_str.split(':')).astype(float))
    delta = -(dec[0]+dec[1]/60.+dec[2]/3600.)
### Reference parameters ###
    tbegin= float(p.get('Reference Parameters','tbegin'))
    talert= float(p.get('Reference Parameters','talert'))
    t0par = float(p.get('Reference Parameters','t0par'))
    t0base= float(p.get('Reference Parameters','t0base'))
    mjd   = p.get('Reference Parameters','mjd')=='True'
### Fitting parameters ###
    t0  = float(p.get('Fitting Parameters','t0'))
    te  = float(p.get('Fitting Parameters','te'))
    u0  = float(p.get('Fitting Parameters','u0'))
    rho = float(p.get('Fitting Parameters','rho'))
    pien= float(p.get('Fitting Parameters','pien'))
    piee= float(p.get('Fitting Parameters','piee'))
### fs,fb of each observatory ###
    I_minus_L = p.get('Flux Parameters','i_minus_l')
    if I_minus_L != '':
        I_minus_L = np.array(I_minus_L.split()).astype(float)
    use_color = p.get('Flux Parameters','use_color')=='True'
    parmflx = []
    errfacs = []
    nameflx = []
    for iob in range(nobs):
        fsfb = p.get('Flux Parameters','(fs,fb)_'+obsnames[iob])
        fsfb = np.array(fsfb.split()).astype(float)
        errfac = p.get('Error Rescalings','errfac_'+obsnames[iob])
        parmflx.append(fsfb)
        errfacs.append(float(errfac))
        nameflx.append(['fs%d'%iob,'fb%d'%iob])
    parmref = np.array([tbegin,talert,t0par,t0base,mjd])
    parmfit = np.array([t0+t0base,u0,te,rho,pien,piee])
    parmflx = np.array(parmflx)
    errfacs = np.array(errfacs)
    namefit = np.array(['t0','u0','te','rho','pien','piee'])
    nameflx = np.array(nameflx)
### print the event information for check ###
    print '####### Event: %s #######'%ename
    print 'RA (J2000): %s'%ra_str,alpha
    print 'Dec (J2000): %s'%dec_str,delta
    print 'nobs = ',nobs
    print 'Data files:'
    print obsnames
    print SPITZ
    print bands
    print datfiles
    print 'Limb-darkening parameters (i,v,h,l): ',gamma_all
    print 'Reference parameters: ',tbegin,talert,t0par,t0base,mjd
    print 'Fitting parameters: '
    print namefit
    print t0+t0base,u0,te,rho,pien,piee
    print 'Flux parameters: '
    print nameflx.flatten()
    print parmflx.flatten()
    print 'use color information? ',use_color,I_minus_L
    print 'Error factors: ',errfacs
    print '##########################'
    return alpha,delta,datfiles,oblats,oblongs,obgammas,SPITZ,parmref,parmfit,parmflx,errfacs,namefit,nameflx,use_color,I_minus_L

if __name__ == '__main__':
    getEventInfo('ob140124')
    getEventInfo('ob140939')
