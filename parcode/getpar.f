c      subroutine getmag_plfs(
      subroutine geta(qn,qe,hjd,alpha,delta,t0)
cf2py intent(in) :: hjd
cf2py intent(in) :: alpha
cf2py intent(in) :: delta
cf2py intent(in) :: t0
cf2py intent(out) :: qn
cf2py intent(out) :: qe
      implicit real*8 (a-h,o-z)
      real*8 sun(3),xpos(3),ypos(3),rad(3),north(3),east(3)
      real*8 spring(3),summer(3)
      data spring/1.,0.,0./
      data summer/0.,0.9175,0.3978/
      data pi/3.14159265/
! warning: data perielio alquanto remota..... (vernal in subroutine(s) get*) [?!]
      common/vearth/vne,vee
c-----------
      ecc = 0.0167
      vernal = 2719.55d0
      offset = 75
      peri   = vernal - offset
      phi = (1 - offset/365.25)*2*pi
      call getpsi(psi,phi,ecc)
      costh = (cos(psi) - ecc)/(1-ecc*cos(psi))
      sinth = -sqrt(1-costh**2)
      do 3 i = 1,3
         xpos(i) = spring(i)*costh + summer(i)*sinth
         ypos(i) =-spring(i)*sinth + summer(i)*costh
 3    continue
 4    format(3f10.4)
      north(1) = 0
      north(2) = 0
      north(3) = 1
      radian = 180/3.14159265
      rad(1) = cos(alpha/radian)*cos(delta/radian)
      rad(2) = sin(alpha/radian)*cos(delta/radian)
      rad(3) = sin(delta/radian)
      call cross(east,north,rad)
      call dot(e2,east,east)
      do 5 i=1,3
         east(i) = east(i)/sqrt(e2)
 5    continue
      call cross(north,rad,east)
 6    format(3f7.3)
      phi   = (t0+1 - peri)/365.25*2*pi
      call getpsi(psi,phi,ecc)
      qn2 = 0
      qe2 = 0
      do 10 i=1,3
         sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
         qn2 = qn2 + sun(i)*north(i)
         qe2 = qe2 + sun(i)*east(i)
 10   continue
      phi   = (t0-1 - peri)/365.25*2*pi
      call getpsi(psi,phi,ecc)
      qn1 = 0
      qe1 = 0
      do 20 i=1,3
         sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
         qn1 = qn1 + sun(i)*north(i)
         qe1 = qe1 + sun(i)*east(i)
 20   continue
      phi   = (t0 - peri)/365.25*2*pi
      call getpsi(psi,phi,ecc)
      qn0 = 0
      qe0 = 0
      do 30 i=1,3
         sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
         qn0 = qn0 + sun(i)*north(i)
         qe0 = qe0 + sun(i)*east(i)
 30   continue
      vn0 = (qn2-qn1)/2
      ve0 = (qe2-qe1)/2
      factor = 365.25*4.74
      vne = -vn0*factor
      vee = -ve0*factor
c	read(5,*)xyz
c  0.0389600164 -0.156597536  0.485358298 -28.9469337
      t = hjd
      phi   = (t - peri)/365.25*2*pi
      call getpsi(psi,phi,ecc)
      qn = -qn0 - vn0*(t-t0)
      qe = -qe0 - ve0*(t-t0)
      do 40 i=1,3
         sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
         qn = qn + sun(i)*north(i)
         qe = qe + sun(i)*east(i)
 40   continue
 11   format(i6,2f9.5)
 100  continue
      return
      end

      subroutine gett(qn,qe,hjd,qlat,qlong,alpha,delta)
cf2py intent(in) :: hjd
cf2py intent(in) :: qlat
cf2py intent(in) :: qlong
cf2py intent(in) :: alpha
cf2py intent(in) :: delta
cf2py intent(out) :: qn
cf2py intent(out) :: qe
      implicit real*8 (a-h,o-z)
      real*8 sun(3),xpos(3),ypos(3),rad(3),north(3),east(3)
      real*8 spring(3),summer(3)
      data pi/3.14159265d0/
c      write(6,*) hjd
      radian = 180/pi
      rearth = 20000/pi/1.5e8 !/100   (remove "/100" to implement)
      vernal = 2719.0d0 + 2.4d0/360.d0
      north(1) = 0
      north(2) = 0
      north(3) = 1
      rad(1) = cos(alpha/radian)*cos(delta/radian)
      rad(2) = sin(alpha/radian)*cos(delta/radian)
      rad(3) = sin(delta/radian)
      call cross(east,north,rad)
      call dot(e2,east,east)
      do 5 i=1,3
         east(i) = east(i)/sqrt(e2)
 5    continue
      call cross(north,rad,east)
c
      qn = 0
      qe = 0
      qr = 0
      dsyn = 365.25d0 - 0.008d0
      phase = (hjd-vernal)*(1d0 + 1d0/dsyn) + qlong/360
c        phase = (hjd-vernal)*366.25/365.25 + qlong/360
        sun(1) = -cos(phase*2*pi)*cos(qlat/radian)
        sun(2) = -sin(phase*2*pi)*cos(qlat/radian)
        sun(3) = -sin(qlat/radian)
      do 30 i=1,3
         qn = qn + sun(i)*north(i)*rearth
         qe = qe + sun(i)*east(i)*rearth
         qr = qr + sun(i)*rad(i)*rearth
 30   continue
      return
      end

      subroutine gets(qn,qe,qr,hjd,
     *   hjdspitz,raspitz,decspitz,disspitz,nspitzer,alpha,delta)
cf2py intent(in) :: hjd
cf2py intent(in) :: hjdspitz
cf2py intent(in) :: raspitz
cf2py intent(in) :: decspitz
cf2py intent(in) :: disspitz
!! nspitzer is optional
cf2py intent(in) :: nspitzer
cf2py intent(in) :: alpha
cf2py intent(in) :: delta
cf2py intent(out) :: qn
cf2py intent(out) :: qe
cf2py intent(out) :: qr
      implicit real*8 (a-h,o-z)
      dimension raspitz(nspitzer),decspitz(nspitzer),disspitz(nspitzer)
      dimension hjdspitz(nspitzer)
      real*8 sun(3),xpos(3),ypos(3),rad(3),north(3),east(3)
      real*8 spring(3),summer(3)
      data pi/3.14159265d0/
      if(hjd.gt.hjdspitz(nspitzer)+1.or.
     *    hjd.lt.hjdspitz(1)-1)then
          qn = 0
          qe = 0
          qr = 0
c          write(6,*)'date out of range: ',hdj,nspitzer,hjdspitz(1)-1
          return
      endif
c   find interpolated ra,dec,dis
      imin1 = 0
      do iloop = 0,1
         dtmin = 1e8
         do i=1,nspitzer
            dt = abs(hjd-hjdspitz(i))
c	      write(6,*)i,imin,dt,dtmin,imin1
          if(dt.lt.dtmin.and.i.ne.imin1)then
    	 imin = i
    	 dtmin = dt
          endif
      enddo
      if(iloop.eq.0)then
          imin1 = imin
      else
          imin2 = imin
      endif
      enddo
      x = (hjd - hjdspitz(imin1))/(hjdspitz(imin2) - hjdspitz(imin1))
      ras = x*raspitz(imin2)   +  (1d0-x)*raspitz(imin1)
      decs = x*decspitz(imin2)   +  (1d0-x)*decspitz(imin1)
      diss = x*disspitz(imin2)   +  (1d0-x)*disspitz(imin1)
c	write(6,*)imin1,imin2,ras,decs,diss
      radian = 180/pi
      north(1) = 0
      north(2) = 0
      north(3) = 1
      rad(1) = cos(alpha/radian)*cos(delta/radian)
      rad(2) = sin(alpha/radian)*cos(delta/radian)
      rad(3) = sin(delta/radian)
      call cross(east,north,rad)
      call dot(e2,east,east)
      do 5 i=1,3
         east(i) = east(i)/sqrt(e2)
  5   continue
      call cross(north,rad,east)
      qn = 0
      qe = 0
      qr = 0
      phase = (hjd-vernal)*366.25/365.25 + qlong/360
c        write(6,*)phase
      sun(1) = -cos(ras/radian)*cos(decs/radian)
      sun(2) = -sin(ras/radian)*cos(decs/radian)
      sun(3) = -sin(decs/radian)
      do 30 i=1,3
         qn = qn + sun(i)*north(i)*diss
         qe = qe + sun(i)*east(i)*diss
         qr = qr + sun(i)*rad(i)*diss
  30  continue
      return
      end

      subroutine getpsi(psi,phi,ecc)
      implicit real*8 (a-h,o-z)
      pi = 3.14159265d0
      psi= phi
      do 10 i=1,4
         fun = psi - ecc*sin(psi)
         dif = phi - fun
         der = 1 - ecc*cos(psi)
         psi = psi + dif/der
 10   continue
      return
      end

      subroutine cross(c,a,b)
      implicit real*8 (a-h,o-z)
      dimension a(3),b(3),c(3)
      c(1) = a(2)*b(3) - b(2)*a(3)
      c(2) = a(3)*b(1) - b(3)*a(1)
      c(3) = a(1)*b(2) - b(1)*a(2)
      return
      end

      subroutine dot(c,a,b)
      implicit real*8 (a-h,o-z)
      dimension a(3),b(3)
      c = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end

