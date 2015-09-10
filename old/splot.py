import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator,FormatStrFormatter
import matplotlib.cm as cm
import getdat
import getparm
import genlc
import sfitpy

def getchi2(nobs,flux,fmod,ferr):
    chi2,dof = [],[]
    for i in range(nobs):
        dof.append(len(flux[i]))
        chi2.append(np.sum((flux[i]-fmod[i])**2/ferr[i]**2))
    return np.array(chi2),np.array(dof)

def main():
    eventname = 'ob150448'
#    eventname = raw_input('Please give the event name (eg, ob140124): ')
    alpha,delta,datfiles,oblats,oblongs,obgammas,isspitz,parmref,parmfit,parmflx,errfacs,namefit,nameflx,use_color,I_minus_L = getparm.getEventInfo(eventname)
    data = getdat.getalldat(datfiles,errfacs,isspitz,parmref)
    tbegin,talert,t0par,t0base,mjd = parmref
    parms_names = [parmfit,namefit,parmflx,nameflx,parmref]
    parmfree,namefree,parmfix,namefix = sfitpy.initialization(parms_names,isspitz)
    ## generating fake data ##
    q,s,theta = 1e-6,10,0.
    parmbfit = np.hstack([parmfit,[theta,q,s]])
    parmbflx = parmflx
#    data_fake = []
#    for iob in range(len(data)):
#        date,flux,ferr,mag,merr = data[iob]
#        ibinlc,xcau,ycau = genlc.getbinlc(alpha,delta,oblats[iob],oblongs[iob],date,parmbfit,parmbflx[iob],obgammas[iob],t0par,isspitz[iob],True)
#        imag = 18-2.5*np.log10(ibinlc[:,2])
#        ## save the fake data sets ##
#        np.savetxt('fake-data-%d.dat'%(iob+1),np.vstack([date,imag,merr,np.ones_like(date),np.ones_like(date)]).T,fmt='%f')
#        data_fake.append([date,ibinlc[:,2],ferr,imag,merr])
#    data = data_fake
#   ## change (pien,piee) ##
#    parmfree[namefree=='pien'] = 0.2
#    parmfree[namefree=='piee'] = 1.
#    print parmfree
    ### Downhill simplex method ###
#    parmfit,parmflx,chi2s,dof = sfitpy.downhill(sfitpy.lnprob_func,namefree,parmfree,namefix,parmfix,data,[alpha,delta,oblats,oblongs,obgammas],parmref,isspitz,use_color=use_color,I_minus_L=I_minus_L)
    ### MCMC to minimize the chisq ###
#    parmfit,parmflx,chi2s,dof = sfitpy.mcmc(sfitpy.lnprob_func,eventname,namefree,parmfree,namefix,parmfix,data,[alpha,delta,oblats,oblongs,obgammas],parmref,isspitz,use_color=use_color,I_minus_L=I_minus_L)
#### use getchi2 instead of lnprob_func: TBD ###
    ### Downhill simplex method ###
    parmfit,parmflx,chi2s,dof = sfitpy.downhill(sfitpy.getchi2,namefree,parmfree,namefix,parmfix,data,[alpha,delta,oblats,oblongs,obgammas],parmref,isspitz,use_color=use_color,I_minus_L=I_minus_L)
    ### MCMC to minimize the chisq ###
#    parmfit,parmflx,chi2s,dof = sfitpy.mcmc(sfitpy.getchi2,eventname,namefree,parmfree,namefix,parmfix,data,[alpha,delta,oblats,oblongs,obgammas],parmref,isspitz,use_color=use_color,I_minus_L=I_minus_L)
#### Find the model according to the given parameters ###
    fmod,mmod = [],[] ## model in flux and magnitude
    trajmod = []
    fmod_fin,mmod_fin = [],[]
    trajmod_fin = []
#    parmfit,parmflx = sfitpy.getparmfit(parmfree,parmfix,namefree,namefix)
    te = parmfit[namefit=='te']
    npts = (7300-6500)*10000./te
    npts = min([npts,100000])
    tmod = np.linspace(6500,7300,int(npts))
    for iob in range(len(data)):
        ## find model ##
        date,flux,ferr = data[iob][:3]
        itraj,ifmod = genlc.getslc(alpha,delta,oblats[iob],oblongs[iob],date,parmfit,parmflx[iob],obgammas[iob],t0par,isspitz[iob],False,True)
        immod = 18.-2.5*np.log10(ifmod)
        fmod.append(ifmod)
        mmod.append(immod)
        trajmod.append(itraj)
        ## find a finer light curve ##
        date = tmod
        itraj,ifmod = genlc.getslc(alpha,delta,oblats[iob],oblongs[iob],date,parmfit,parmflx[iob],obgammas[iob],t0par,isspitz[iob],False,True)
        immod = 18.-2.5*np.log10(ifmod)
        fmod_fin.append(ifmod)
        mmod_fin.append(immod)
        trajmod_fin.append(itraj)
#    fmod = np.array(fmod)
#    mmod = np.array(mmod)
    trajmod = np.array(trajmod)
    trajmod_fin = np.array(trajmod_fin)
    print 'Total chisq = ',chi2s[0]
    print ['%.2f/%d'%(chi2s[i+1],dof[i]) for i in range(len(dof))]
### start plotting ###
    fig = plt.figure(figsize=(8.5,11))
    ax = plt.subplot(211)
    ax_traj = plt.subplot(212)
    ax.set_title(eventname)
#    colorstr = ['k','r','g','m']
    colorstr = cm.rainbow(np.linspace(0,1,len(data)))
    ground_mod = True
    for iob in range(len(data)):
        ## plot data ##
        date,flux,ferr = data[iob][:3]
        imag,imerr = getdat.convert_to_ogle_dat(flux,ferr,parmflx[iob],parmflx[0],use_mag=True)
        ax.errorbar(date,imag,yerr=imerr,marker='o',linestyle='none',color=colorstr[iob],markeredgecolor=colorstr[iob],markerfacecolor='none')
        ## plot model ##
        imod = getdat.convert_to_ogle_mod(fmod[iob],parmflx[iob],parmflx[0],use_mag=True)
        imod_fin = getdat.convert_to_ogle_mod(fmod_fin[iob],parmflx[iob],parmflx[0],use_mag=True)
        np.savetxt('output/%s-%d.dat'%(eventname,iob+1),np.vstack([date,trajmod[iob].T,imag,imerr,imod]).T,fmt='%f')
        np.savetxt('output/%s-%d-fin.dat'%(eventname,iob+1),np.vstack([tmod,trajmod_fin[iob].T,imod_fin]).T,fmt='%f')
        if isspitz[iob]==False and ground_mod==False:
            continue
        ax.plot(tmod,imod_fin,color=colorstr[iob])
        ax_traj.plot(trajmod[iob][:,0],trajmod[iob][:,1],marker='o',color=colorstr[iob])
        ax_traj.plot(trajmod_fin[iob][:,0],trajmod_fin[iob][:,1],color=colorstr[iob])
        ground_mod = False
    ax.invert_yaxis()
    ax.set_xlabel('HJD')
    ax.set_ylabel('I (OGLE)')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax_traj.plot(0,0,marker='o',markersize=4,color='k')
    ax_traj.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_traj.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax_traj.axis('equal')
    ax_traj.set_xlabel('x')
    ax_traj.set_ylabel('y')
    plt.show()
    return

if __name__ == '__main__':
    main()
