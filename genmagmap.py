import numpy as np
import caustic_double_precision as gencau
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import os
import getbmag
import getmag_wei
import time

def mapmaking(q,s,rho,Gamma,maplimits,gridsizes):
    xmin,xmax,ymin,ymax = maplimits
    xstep,ystep = gridsizes
    rho2 = rho**2
    ## lens positions centered on the center of magnification ##
    if s<1:
        offset_magcen = q/(1+q)*s
    else:
        offset_magcen = q/(1+q)/s
    x1 = -offset_magcen
    x2 = s-offset_magcen
#    zcau,zcr = gencau.getCaustic(np.array([1.,q]),np.array([0.,s]),np.array([0.,0.]),1000)
    zcau,zcr = gencau.getCaustic(np.array([1.,q]),np.array([x1,x2]),np.array([0.,0.]),1000)
    xcau = np.real(zcau)
    ycau = np.imag(zcau)
    if s < 1.:
        offset_x = s*(-0.5+1./(1+q))
    else:
        offset_x = s/2.-q/(1+q)/s
#    offset_x = s/2.
    magmap = []
    for xs in np.arange(xmin,xmax,xstep):
        for ys in np.arange(ymin,ymax,ystep):
            dist2 = min((xcau-xs)**2+(ycau-ys)**2)
            xs_corr = xs-offset_x
            if dist2 < 4*rho2:
            ### Jin code ###
                magbps,magnew,errorflag = getbmag.getmag_jin(xs_corr,ys,s,q,rho,Gamma)
                if magnew < 1.:
                    print 'looplinking failed, switch to mapmaking'
                    magnew,errorflag = getmag_wei.getmag_wei(magbps,xs_corr,ys,s,q,rho,Gamma)
            elif dist2 < 16*rho2:
            ### hex approx ###
                magnew = getbmag.taylor_2(xs_corr,ys,q,s,Gamma,rho,4)
            elif dist2 < 100*rho2:
            ### quad approx ###
                magnew = getbmag.taylor_2(xs_corr,ys,q,s,Gamma,rho,2)
            else:
            ### mon approx ###
                magnew = getbmag.taylor_2(xs_corr,ys,q,s,Gamma,rho,0)
            magmap.append([xs,ys,magnew])
    return np.array(magmap),xcau,ycau

def pltmap(ax,x,y,z,rho,maplimits,gridsizes):
    xmin,xmax,ymin,ymax = maplimits
    xstep,ystep = gridsizes
    xi = np.arange(xmin,xmax,xstep)
    yi = np.arange(ymin,ymax,xstep)
    zi = griddata(x,y,np.log10(z),xi,yi,interp='linear')
    ax.contourf(xi,yi,zi,len(xi),cmap='gray')
#    ax.contourf(xi,yi,zi,len(xi),cmap='OrRd')
#    cbar = plt.colorbar()
    return

def main():
    q = 0.5
    s = 1.5
    rho = 0.01
    xmin,xmax,ymin,ymax = -0.1,1.0,-0.5,0.5
    xstep,ystep = 0.5*rho,0.5*rho
    
    maplimits = [xmin,xmax,ymin,ymax]
    gridsizes = [xstep,ystep]
    Gamma = 0.
    start = time.clock()
    magmap,xcau,ycau = mapmaking(q,s,rho,Gamma,maplimits,gridsizes)
    print 'mapmaking finished: ',time.clock()-start
    print np.shape(magmap)
    x = magmap[:,0]
    y = magmap[:,1]
    mag = magmap[:,2]
    print max(mag),min(mag)
    
    ax = plt.subplot(111)
    ax.plot(xcau,ycau,linestyle='none',marker='o',markersize=1,color='k')
    pltmap(ax,x,y,mag,rho,maplimits,gridsizes)
    ax.axis('equal')
    plt.show()
    return

if __name__ == '__main__':
    main()
