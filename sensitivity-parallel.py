import numpy as np
#import genmagmap
import getdat
import getparm
import genlc
import sfitpy
import time
import multiprocessing as mp

def get_q_s_prob(q,s,eventname,eventinfo,data,parms_names,isspitz,use_color,I_minus_L,ntheta,outputlc=False,finmod_opt=False,date_mod=[]):
    alpha,delta,oblats,oblongs,obgammas = eventinfo
    parmfit,namefit,parmflx,nameflx,parmref = parms_names
    tbegin,talert,t0par,t0base,mjd = parmref
    thetas = np.linspace(0,360,ntheta+1)
    binlcs,binfull = [],[]
    bincau = []
    use_fsfb = True
    start = time.clock()
    ############################
    ### fake data generation ###
    ############################
    print '----- Generating binary lcs... ------'
    for i in range(ntheta):
        itheta = thetas[i]
        parmbfit = np.hstack([parmfit,[itheta,q,s]])
        parmbflx = parmflx
        iobbinlc,iobfull = [],[]
        for iob in range(len(data)):
            date,flux,ferr,mag,merr = data[iob]
            ibinlc,xcau,ycau = genlc.getbinlc(alpha,delta,oblats[iob],oblongs[iob],date,parmbfit,parmbflx[iob],obgammas[iob],t0par,isspitz[iob],use_fsfb,caustic=bincau)
            bincau = [xcau,ycau]
            iobbinlc.append([date,ibinlc[:,2],ferr])
            iobfull.append([date,ibinlc[:,0],ibinlc[:,1],ibinlc[:,2],ferr])
            ### find the fine input model ###
            if finmod_opt == True:
                ibinlc_mod,xcau,ycau = genlc.getbinlc(alpha,delta,oblats[iob],oblongs[iob],date_mod,parmbfit,parmbflx[iob],obgammas[iob],t0par,isspitz[iob],use_fsfb,caustic=bincau)
                np.savetxt('output/models/finmod-%s-%.5f-%.2f-%03d-%d.dat'%(eventname,q,s,itheta,iob+1),np.vstack([date_mod,ibinlc_mod.T]).T,fmt='%f')
        binlcs.append(iobbinlc)
        binfull.append(iobfull)
    np.savetxt('output/caustics/caustic-%s-%.5f-%.2f.dat'%(eventname,q,s),np.array(bincau).T,fmt='%f')
    print '------ Binary lcs generated: ',time.clock()-start,' s ------'
    ##################################
    ### employ sfitpy to fake data ###
    ##################################
    if talert > max(data[0][0]): ## obj event
        print '--- This is an obj event ---'
        trials = sfitpy.massfit_obj(q,s,thetas,binlcs,binfull,eventinfo,parms_names,isspitz,use_color,I_minus_L,outputlc)
    else: ## subj event
        print '--- This is a subj event ---'
        trials = sfitpy.massfit_sub(q,s,thetas,binlcs,binfull,eventinfo,parms_names,isspitz,use_color,I_minus_L,outputlc)
    return trials


def main(eventname):
    alpha,delta,datfiles,oblats,oblongs,obgammas,isspitz,parmref,parmfit,parmflx,errfacs,namefit,nameflx,use_color,I_minus_L = getparm.getEventInfo(eventname)
    data = getdat.getalldat(datfiles,errfacs,isspitz,parmref)
    tbegin,talert,t0par,t0base,mjd = parmref
### use more intensive time sequence ###
    date_mod = np.linspace(6700,7000,2000)
##### choose s and q for sensitivity mapping ###
#    ## for ob0124 ##
#    ntheta = 300
#    qList = np.logspace(np.log10(1e-5),np.log10(1e-2),20)
#    sList = np.logspace(-0.6,0.6,30)
#    nprocesses = 25
#    ## for ob0939 ##
#    ntheta = 300
#    qList = np.logspace(np.log10(1e-4),np.log10(1e-2),20)
#    sList = np.logspace(-0.4,0.4,100)
#    nprocesses = 6
#    ### pool.apply_async to parallelize ###
#    pool = mp.Pool(nprocesses)
#    for q in qList:
#        for s in sList:
#            pool.apply_async(get_q_s_prob,args=(q,s,eventname,[alpha,delta,oblats,oblongs,obgammas],data,[parmfit,namefit,parmflx,nameflx,parmref],isspitz,use_color,I_minus_L,ntheta),callback=mycallback)
#    pool.close()
#    pool.join()
#################################################
#### for test purpose ###
    plot_opt = True
    finmod_opt = True
    ntheta = 5
    qList = [7e-4]
    sList = [0.94]
#######################
    for q in qList:
        for s in sList:
            trials = get_q_s_prob(q,s,eventname,[alpha,delta,oblats,oblongs,obgammas],data,[parmfit,namefit,parmflx,nameflx,parmref],isspitz,use_color,I_minus_L,ntheta,outputlc=True,finmod_opt=True,date_mod=date_mod)
            mycallback(trials)
    return

def mycallback(trials):
    global outputs
    if outputs == []:
        outputs = trials
    else:
        outputs = np.vstack([outputs,trials])

eventname = 'ob150448'
outputs = []
main(eventname)
np.savetxt('output/grid-%s.dat'%eventname,outputs,fmt='%e')
