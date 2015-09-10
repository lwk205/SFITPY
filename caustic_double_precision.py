import numpy as np
import scipy as sp
#import matplotlib.pyplot as plt
import time

def getCaustic(m, x, y, N):
    nlens = len(x)
    z = x+y*1j
    totalMass = sum(m)
#    totalMass = max(m)
    m /= totalMass
    
    f0 = np.zeros(2*nlens+1)*1j
    g = np.zeros([nlens,2*nlens])*1j
    f = np.zeros([nlens,2*nlens])*1j
    dphi = 2.0*np.pi/N
    zsList = []
    zcrList = []
    for parts in range(N):
    	phi = parts*dphi
    	f0[0] = z[0]**2
    	f0[1] = -2.0*z[0]
    	f0[2] = 1.0
    	k = 1
    	for lens in range(1,nlens):
    	    k += 2
    	    f0[k+1] = f0[k-1]
    	    f0[k] = f0[k-2]-2.0*f0[k-1]*z[lens]
    	    for j in range(k-2,0,-1):
    		f0[j+1] = f0[j-1]-2.0*f0[j]*z[lens]+f0[j+1]*z[lens]**2
    	    f0[1] = -2.0*f0[0]*z[lens]+f0[1]*z[lens]**2
    	    f0[0] = f0[0]*z[lens]**2
    	for lens in range(nlens):
    	    g[lens,2*nlens-1] = f0[2*nlens]
    	    for j in range(2*nlens,1,-1):
    		g[lens,j-2] = g[lens,j-1]*z[lens]+f0[j-1]
    	    f[lens,2*nlens-2] = g[lens,2*nlens-1]
    	    for j in range(2*nlens-1,1,-1):
    		f[lens,j-2] = f[lens,j-1]*z[lens]+g[lens,j-1]
    
        c = np.zeros(2*nlens+1)*1j
    	eiphi = complex(np.cos(phi), np.sin(phi))
    	c[2*nlens] = f0[2*nlens]*eiphi
    	c[2*nlens-1] = f0[2*nlens-1]*eiphi
    	for order in range(2*nlens-1,0,-1):
    	    c[order-1] = f0[order-1]*eiphi
    	    secondTerm = 0.0
    	    for j in range(nlens):
    		secondTerm += m[j]*f[j,order-1]
    	    c[order-1] -= secondTerm
    
    	orders = 2*nlens+1
        coeff = np.zeros(orders)*1j
    	for order in range(orders):
    	    coeff[order] = c[orders-order-1]
    
    	zcr = np.roots(coeff)
        zs = np.zeros_like(zcr)*1j
    	for eachRoot in range(len(zcr)):
    	    zs[eachRoot] = zcr[eachRoot]
    	    for lens in range(nlens):
    		zs[eachRoot] -= m[lens]/(np.conj(zcr[eachRoot])-np.conj(z[lens]))
    
        zsList.extend(zs)
        zcrList.extend(zcr)
    return np.array(zsList),np.array(zcrList)

def main():
    q = 7e-4
    s = 0.94
    ## center of magnification ##
    if s<1:
        offset_magcen = q/(1+q)*s
    else:
        offset_magcen = q/(1+q)/s
    x1 = -offset_magcen
    x2 = s-offset_magcen
    m = np.array([1.0,q])
    x = np.array([x1,x2])
    y = np.array([0.,0.])
    zcau,zcr = getCaustic(m,x,y,10000)
    xcau = np.real(zcau)
    ycau = np.imag(zcau)
    np.savetxt('caustic.dat',np.vstack([xcau,ycau]).T,fmt='%f')
    return

if __name__ == '__main__':
    main()
