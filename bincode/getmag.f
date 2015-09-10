      Include 'getbinp.f'
      Include 'jinsubs.f'
      subroutine taylor_2(qmag, x, y, q, d, gamma,
     + rho,ntaylor)
cf2py intent(in) :: x
cf2py intent(in) :: y
cf2py intent(in) :: q
cf2py intent(in) :: d
cf2py intent(in) :: gamma
cf2py intent(in) :: rho
cf2py intent(in) :: ntaylor
cf2py intent(out) :: qmag
      implicit real*8 (a-h,o-z)
      real*8 q, d
      integer nsol
      integer MAX_SOL
      parameter(MAX_SOL = 5)
      real*8 xic(MAX_SOL), yic(MAX_SOL), am(MAX_SOL)
ccc      call
ccc     + getimt(xic, yic, am, 
ccc     + nsol, da, db, phi, ma, mb, x, y, qmag)
      call getbinp(qmag,am,xic,yic,x,y,q,d,nsol)
      if(ntaylor.eq.0)return
      qmag0 = qmag
      call getquad(qmag1,x,y,q,d,rho,0)
      qmag1 = qmag1 - qmag0
      if(ntaylor.eq.2)then
         qmag = qmag0 + qmag1/2.d0*(1.d0 - gamma/5.d0)
         return
      endif
      if(ntaylor.eq.4)then
         call getquad(qmag1p,x,y,q,d,rho,1)
        qmag1p = qmag1p - qmag0
         call getquad(qmaghalf,x,y,q,d,rho/2,0)
         qmaghalf = qmaghalf - qmag0
         a = (16*qmaghalf - qmag1)/3
         b = (qmag1 + qmag1p)/2 - a
         qmag = qmag0 + a/2.d0*(1.d0 - gamma/5.d0)
     *        + b/3.d0*(1.d0 - gamma*11.d0/35.d0)
         return
      endif
      write(6,*)ntaylor,' bad ntaylor'
      stop
      end


      subroutine getquad(qmag1,x0,y0,q,b,rho,nrot)
      implicit real*8 (a-h,o-z)
      dimension npos(5)
      data npos/1,0,-1,0,1/
      real*8 q, b 
      integer nsol
      integer MAX_SOL
      parameter(MAX_SOL = 5)
      real*8 xic(MAX_SOL), yic(MAX_SOL), am(MAX_SOL)
      qmag1 = 0
      if(nrot.eq.0)then
         do 10 i=1,4 
            x = x0 + npos(i)*rho
            y = y0 + npos(i+1)*rho
c            call getpt(qmag,x,y,da,db,ma,mb,phi)
ccc       call
ccc     + getimt(xic, yic, am,
ccc     + nsol, da, db, phi, ma, mb, x, y, qmag)
       call getbinp(qmag,am,xic,yic,x,y,q,b,nsol)

            qmag1 = qmag1 + qmag
 10      continue
         qmag1 = qmag1/4d0
         return
      endif
      if(nrot.eq.1)then
         fac = sqrt(0.5d0)
         do 20 i=-1,1,2
            do 20 j=-1,1,2
               x = x0 + i*rho*fac
               y = y0 + j*rho*fac
c               call getpt(qmag,x,y,da,db,ma,mb,phi)
cc       call
cc     + getimt(xic, yic, am,
cc     + nsol, da, db, phi, ma, mb, x, y, qmag)
       call getbinp(qmag,am,xic,yic,x,y,q,b,nsol)
               qmag1 = qmag1 + qmag
 20      continue
         qmag1 = qmag1/4
         return
      endif
      write(6,*)nrot,' bad nrot'
      stop
      end


