import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import genmagmap
import matplotlib.pyplot as plt
import getdat
import getparm
import genlc
import sfit
import time

def getchi2(nobs,flux,fmod,ferr):
    chi2,dof = [],[]
    for i in range(nobs):
        dof.append(len(flux[i]))
        chi2.append(np.sum((flux[i]-fmod[i])**2/ferr[i]**2))
    return np.array(chi2),np.array(dof)


def main():
    eventname = 'ob140939'
### find event information, input data, program setup ###
    alpha,delta,datfiles,oblat,oblong,obgamma = getparm.getEventInfo(eventname)
    parmref,parmfit,parmflx,efacs = getparm.getFitSetting(eventname)
    t0par = parmref[0]
    t0base= parmref[1]
    HMJD  = parmref[2]
    nobs  = len(datfiles)
    parmfit[0] += t0base
    print 'Total number of obs: ',nobs
    print 'Parameters (ref): ',parmref
    print 'Parameters (sfit): ',parmfit
    print 'Parameters (flux): ',parmflx
    print 'Errfacs: ',efacs
#    isspit = np.zeros(nobs).astype(bool)
    isspit = np.zeros(nobs-1).astype(bool)
    isspit = np.hstack([isspit,True])
    linestrs = ['-','--']
    date,flux,ferr,magr,merr = getdat.getdat1(nobs,datfiles,efacs,isspit,HMJD)
#############################################################
### use more intensive time sequence ###
    date_mod = np.linspace(6800,6900,1000)
#### choose s and q for sensitivity mapping ###
#    plot_opt = False
#    finmod_opt = False
#    ntheta = 300
#    qList = np.logspace(np.log10(1e-4),np.log10(1e-2),50)
#    sList = np.logspace(np.log10(0.1),np.log10(10),50)
#    ### five parallel programs ###
#    ncore = 6
#    each = 50/ncore
##    qList = qList[:each]
##    qList = qList[each:2*each]
##    qList = qList[2*each:3*each]
##    qList = qList[3*each:4*each]
##    qList = qList[4*each:5*each]
#    qList = qList[5*each:]
#################################################
### for test purpose ###
    plot_opt = True
    finmod_opt = True
    ntheta = 6
    qList = [1e-2]
    sList = [1.1]
########################
    for q in qList:
        trials = []
        for s in sList:
            theta = np.linspace(0,360,ntheta+1)
            xtrajs,ytrajs,binlcs = [],[],[]
            binlcs_mod = []
            use_fsfb = True
            start = time.clock()
            ############################
            ### fake data generation ###
            ############################
            for i in range(ntheta):
                itheta = theta[i]
                parmbfit = np.hstack([parmfit,[itheta,q,s]])
                parmbflx = parmflx
                iobxtraj,iobytraj,iobbinlc = [],[],[]
                iobbinlc_mod = []
                for iob in range(nobs):
                    ibinlc,xcau,ycau = genlc.getbinlc(alpha,delta,oblat[iob],oblong[iob],date[iob],parmbfit,parmbflx[iob],obgamma[iob],t0par,isspit[iob],use_fsfb)
                    iobxtraj.append(ibinlc[:,0])
                    iobytraj.append(ibinlc[:,1])
                    iobbinlc.append(ibinlc[:,2])
                    ### find the fine input model ###
                    if finmod_opt == True:
                        ibinlc_mod,xcau,ycau = genlc.getbinlc(alpha,delta,oblat[iob],oblong[iob],date_mod,parmbfit,parmbflx[iob],obgamma[iob],t0par,isspit[iob],use_fsfb)
                        iobbinlc_mod.append(ibinlc_mod[:,2])
                xtrajs.append(iobxtraj)
                ytrajs.append(iobytraj)
                binlcs.append(iobbinlc)
                binlcs_mod.append(iobbinlc_mod)
            xtrajs = np.array(xtrajs)
            ytrajs = np.array(ytrajs)
            binlcs = np.array(binlcs)
            print 'binary lcs generated: ',time.clock()-start
########################################################
            ################################
            ### employ sfit to fake data ###
            ################################
########################################################
            sinlcs,chi2all = [],[]
            parmfit0 = parmfit
            parmflx0 = parmflx
            for i in range(ntheta):
                itheta = theta[i]
########################################
############# SFITPY MODULE ############
########################################
### if sfit.getchi2 is used, the plotted
### light curve will not be right, beca-
### -use fs, fb are still the input vals
########################################
######## Downhill minimization only ###
                iparmfit,iparmflx,chi2 = sfit.downhill(sfit.lnprob_func,parmfit,parmflx,[date,binlcs[i],ferr],nobs,[alpha,delta,oblat,oblong,obgamma],t0par,isspit)
#                iparmfit,iparmflx,chi2 = sfit.downhill(sfit.getchi2,parmfit,parmflx,[date,binlcs[i],ferr],nobs,[alpha,delta,oblat,oblong,obgamma],t0par,isspit)
######## MCMC minimization only ###
#                iparmfit,iparmflx,chi2 = sfit.mcmcfit(sfit.getchi2,parmfit,parmflx,[date,binlcs[i],ferr],nobs,[alpha,delta,oblat,oblong,obgamma],t0par,isspit)
######## Downhill + MCMC ####
#                iparmfit,iparmflx,chi2 = sfit.downhill(sfit.getchi2,parmfit,parmflx,[date,binlcs[i],ferr],nobs,[alpha,delta,oblat,oblong,obgamma],t0par,isspit)
#                if chi2[0] > 200:
#                    iparmfit,iparmflx,chi2 = sfit.mcmcfit(sfit.getchi2,iparmfit,iparmflx,[date,binlcs[i],ferr],nobs,[alpha,delta,oblat,oblong,obgamma],t0par,isspit)
########################################
                if chi2[0] < 100:
                    parmfit0 = iparmfit
                    parmflx0 = iparmflx
                ### find best-fit model ###
                isinlc = []
                for iob in range(nobs):
                    itraj,ifmod = genlc.getslc(alpha,delta,oblat[iob],oblong[iob],date_mod,iparmfit,iparmflx[iob],obgamma[iob],t0par,isspit[iob],False,True)
                    isinlc.append(ifmod)
                ###########################
                temp = np.hstack([q,s,theta[i],chi2.reshape(-1),iparmfit,iparmflx.reshape(-1)])
                trials.append(temp)
                sinlcs.append(isinlc)
                chi2all.append(chi2[0])
                print '##########',q,s,theta[i],chi2[0],'############'
            sinlcs = np.array(sinlcs)
############################################################
            ###############################
            ### fake data demonstration ###
            ###############################
            if plot_opt == True:
                color1 = ['#99d8c9','#66c2a4','#41ae76','#238b45','#005824'] #green
                color2 = ['#fdbb84','#fc8d59','#ef6548','#d7301f','#990000'] #red
                color3 = ['#a6bddb','#74a9cf','#3690c0','#0570b0','#034e7b'] #blue
                color4 = ['#bdbdbd','#969696','#737373','#525252','#252525'] #black
                colors = np.array([color1,color2,color3,color4]).reshape(-1)
                fig1 = plt.figure(figsize=(8.5,11))
                nax1 = ntheta/2
                axs = []
                fig2 = plt.figure(figsize=(8.6,8))
                ax0 = fig2.add_subplot(111) # axis to plot the trajectories
                ax0.plot(xcau,ycau,linestyle='none',color='k',marker='o',markersize=1)
                ax0.axis('equal')
                for i in range(ntheta):
                    itheta = theta[i]
                    axi = fig1.add_subplot(nax1,2,i+1)
                    axs.append(axi)
                    for iob in range(nobs):
                        ### plot fake data ###
#                        iflux,iferr = binlcs[i][iob],ferr[iob] #do not convert to OGLE system
                        iflux,iferr = getdat.convert_to_ogle_dat(binlcs[i][iob],ferr[iob],parmflx[iob],parmflx[0]) #convert to OGLE system
                        axi.errorbar(date[iob],iflux,yerr=iferr,color=colors[i],linestyle='none',marker='o')
                        if finmod_opt == True:
#                            ifmod = binlcs_mod[i][iob] #do not convert to OGLE system
                            ifmod = getdat.convert_to_ogle_mod(binlcs_mod[i][iob],parmflx[iob],parmflx[0]) #convert to OGLE system
                            axi.plot(date_mod,ifmod,color=colors[i],linestyle='--')
                        ### plot best-fit smod ###
#                        ifmod = sinlcs[i][iob] #do not convert to OGLE system
                        ifmod = getdat.convert_to_ogle_mod(sinlcs[i][iob],parmflx[iob],parmflx[0]) #convert to OGLE system
                        axi.plot(date_mod,ifmod,color=colors[i])
                        ### plot trajectory ###
                        ax0.plot(xtrajs[i][iob],ytrajs[i][iob],color=colors[i],linestyle=linestrs[iob])
                    axi.set_xlim(6800,6900)
                    x1,x2,y1,y2 = axs[i].axis()
                    xpos = x1+(x2-x1)/10.
                    ypos = y2-(y2-y1)/10.
                    axi.text(xpos,ypos,r'$\theta=%3.1f^\circ, \chi^2=%3.1f$'%(itheta,chi2all[i]),va='center')
                    axi.invert_yaxis()
                plt.show()
        ### save parameters and chi2 ###
        np.savetxt('output/grid-%s-%e.dat'%(eventname,q),np.array(trials),fmt='%f')
    return

if __name__ == '__main__':
    main()
