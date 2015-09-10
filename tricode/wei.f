      include 'roots10.f'
      include '../bincode/getbinp.f'
      include 'getimt2.f'
      include 'triple_lens2.f'
      subroutine getmag_wei(magbfs,magsingle,xs,ys,b,q,rhos,
     + gamma,errorflag)
cf2py intent(out) :: magbfs
cf2py intent(in) :: magsingle
cf2py intent(in) :: xs
cf2py intent(in) :: ys
cf2py intent(in) :: b
cf2py intent(in) :: q
cf2py intent(in) :: rhos
cf2py intent(in) :: gamma
cf2py intent(out) :: errorflag
      real*8 xs,ys,b,q,gamma,rhos
      real*8 magbfs,magsingle
      integer errorflag
      real*8 xorigin,yorigin,loopgridsize
      real*8 fstol !precision for looplinking
      parameter (PI=3.141592653589793238462643d0)
      complex*16 zr(10)

      errorflag = 0
      fstol = 0.001d0
      xorigin = 0.d0
      yorigin = 0.d0
      if (magsingle .gt. 50) magsingle = 50d0
      loopgridsize = (fstol*sqrt(1.5d0*PI*PI*magsingle))**(2d0/3d0)*rhos
c      write(*,*) loopgridsize
      call GETMAG(b,5d0,90d0,q,1e-5,xs,ys,rhos,magbfs,gamma,
     + loopgridsize,xorigin,yorigin,errorflag,zr)
      return 
      END
     
c      subroutine get_t0
c     +           (t0,aoutp,nparm,
c     +           qlati, qlongi, rai, deci, t0par, tbinary)
c      implicit none
c      integer*4 k,ndata,nparm
c      real*8 mapgridsize
cc     loopgridsize
c      real*8 xorigin,yorigin
c      integer*4 errorflag
c      integer error
c      real*8 PI 
c      parameter (PI=3.141592653589793238462643d0)
c      real*8 rhos,magnew
c      real*8 magbps, magbfs
c      real*8 xs,ys,fstol
c      real*8 theta,u0,t,m1,m2
c      real*8 t0,te,tau
c      real*8 aoutp(nparm)
c      real*8 tmin,tmax
c      integer MAX_SOL
c      parameter (MAX_SOL=5)
c      real*8 xic(MAX_SOL), yic(MAX_SOL), am(MAX_SOL)
c      integer*4 nsol
c      integer*4 ntab 
c      real*8 b0,b1,db0,db1
c
c      integer i
c
c      real*8 offset_x, offset_y
cc      real*8 da, db, ma, mb, phi
c      real*8 b, q
c      real*8 db, dep, b_0, ep_0, ep, tbinary, xcm, ycm
c
c      logical magfail
c
c      real*8 qlati, qlongi
c      real*8 qn, qe, piex, piey, qnp, qep, taup, betap
c      real*8 dtau, dbeta
c      real*8 rai, deci, t0par
c
c      real*8 magtemp
c
c      integer nspitzer
c      real*8 t0min,t0max,t0_test,tcc
c      real*8 mag1,mag2
c      integer ii,imax
c      integer nsol1,nsol2,nsol_test
c
cc-----Fixed paramters
c      ep_0=0d0 !reference binary rotation
cc      tbinary = 4716.85d0 !reference time (where the map is calculated)
cc---------------------
c
c      magfail = .false.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      tcc   = aoutp(1)
c      u0    = aoutp(2)
c      te    = aoutp(3)
c      rhos  = aoutp(4)
c      piex  = aoutp(5)
c      piey  = aoutp(6)
c      theta = aoutp(7)
c      db    = aoutp(8)
c      dep   = aoutp(9)
c      b_0   = aoutp(10)
c      q     = aoutp(11)
cc$$$      b     = aoutp(8)
cc$$$      q     = aoutp(9)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc find t0 according to tcc and other parameters ccc
c      t0min = 6820.d0
c      t0max = 6850.d0
c      imax = 1000
c            call getmag_tcc
c     +           (nsol1,aoutp,nparm,t0min,
c     +           qlati, qlongi, rai, deci, t0par,tbinary)
cc            if (nsol1 .ne. 5) then
cc                write(6,*) 't0min too large',nsol1
cc                stop
cc            endif
c            call getmag_tcc
c     +           (nsol2,aoutp,nparm,t0max,
c     +           qlati, qlongi, rai, deci, t0par,tbinary)
c            if (nsol1 .eq. nsol2) then
c                write(6,*) 'nsol1 = nsol2 = ',nsol1
c                stop
c            endif
c      do 909 ii=1,imax
c            t0_test = (t0min+t0max)/2.d0
c            call getmag_tcc
c     +           (nsol_test,aoutp,nparm,t0_test,
c     +           qlati, qlongi, rai, deci, t0par,tbinary)
c            if (nsol_test .eq. nsol1) t0min = t0_test
c            if (nsol_test .eq. nsol2) t0max = t0_test
cc            write(6,*) nsol_test
c            if (abs(t0min-t0max)<1e-4) then
c                goto 910
c            endif
c909     enddo
c910     continue
c        t0 = t0min
c      return
c      END
c
c      subroutine getmag_tcc
c     +           (nsol,aoutp,nparm,t0,
c     +           qlati, qlongi, rai, deci, t0par,tbinary)
c      implicit none
c      integer*4 k,nparm
c      real*8 xorigin,yorigin
c      integer*4 errorflag
c      integer error
c      real*8 PI 
c      parameter (PI=3.141592653589793238462643d0)
c      real*8 rhos,magnew
c      real*8 magbps, magbfs
c      real*8 xs,ys,fstol
c      real*8 theta,u0,t,m1,m2
c      real*8 t0,te,tau,tcc
c      real*8 aoutp(nparm)
c      logical flag,switchflag
c      real*8 magpoint, magfspoint
c      integer MAX_SOL
c      parameter (MAX_SOL=5)
c      real*8 xic(MAX_SOL), yic(MAX_SOL), am(MAX_SOL)
c      integer*4 nsol
c      logical dofs
cc      real*8 tmin,tmax
c      real*8 x,x2 
c      integer*4 ntab 
c      real*8 b0,b1,db0,db1
c      real*8 s_r_max, s_r_max_minus2, rsource
c
c      integer i
c
c      real*8 offset_x, offset_y
c      real*8 b, q
c      real*8 db, dep, b_0, ep_0, ep, tbinary, xcm, ycm
c
c      logical magfail
c
c      real*8 qlati, qlongi
c      real*8 qn, qe, piex, piey, qnp, qep, taup, betap
c      real*8 dtau, dbeta
c      real*8 rai, deci, t0par
c
c      real*8 magtemp
c
c      integer nspitzer
c      real*8 xyz
cc-----Fixed paramters
c      ep_0=0d0 !reference binary rotation
cc---------------------
c      magfail = .false.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      tcc   = aoutp(1) !use tcc as free parameter
c      u0    = aoutp(2)
c      te    = aoutp(3)
c      rhos  = aoutp(4)
c      piex  = aoutp(5)
c      piey  = aoutp(6)
c      theta = aoutp(7)
c      db    = aoutp(8)
c      dep   = aoutp(9)
c      b_0   = aoutp(10)
c      q     = aoutp(11)
c
c      offset_y = 0.d0
c
c      !tcc is defined for CTIO
c            t = tcc
c            tau=(t-t0)/te
c
c            b = b_0+db*(t-tbinary)/365.25
c            ep = ep_0 + dep*(t-tbinary)/365.25
c
c            offset_x = b*(-0.5d0+1d0/(1d0+q))
c
c            call geta(qn, qe, t, rai, deci, t0par)
c            call gett(qnp, qep, t, qlati, qlongi,
c     +           rai, deci)
c
c            qn = qn + qnp
c            qe = qe + qep
c            dtau = piex*qn + piey*qe
c            dbeta = -piex*qe + piey*qn
c            dbeta = -dbeta
c            taup  = tau + dtau
c            betap = u0 + dbeta
c
c            xcm= taup*dcos(theta)+betap*dsin(theta)! -offset_x
c            ycm=-taup*dsin(theta)+betap*dcos(theta)!  -offset_y
c
c            xs = xcm*dcos(ep)-ycm*dsin(ep)-offset_x
c            ys = xcm*dsin(ep)+ycm*dcos(ep)
c
cccc get the point-source magnification at this point
c            call getbinp(magnew,am,xic,yic,xs,ys,q,b,nsol)
cc            write(6,*) tcc,t0,te,nsol,xs,ys,magnew
c        return
c        END
