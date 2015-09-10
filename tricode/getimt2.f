      subroutine getimt(xi,yi,am,nsol,da,db,phi,ma,mb,xs,ys,magtot,zr)
CU    uses lenseqt,lenseqtJ,conj3,dmag,zroots3,laguer3
      double precision PI,RADIAN
      parameter(PI=3.141592653589793238462643d0)
      parameter(RADIAN=1.8d2/PI)
      double precision tol
      double precision da,db,x1,x3,y3,phi
      double precision ma,mb,q1,q2,q3
      double precision xs,ys,magtot,derr,detJ
      double precision xi(10),yi(10),am(10)
      complex*16 z1,z2,z3
      complex*16 a,b,c,r,s,p,q,err
      complex*16 ab,bb,cb,pb,qb
      complex*16 alpha,beta,gamma,delta,eta,zeta,zetab
      complex*16 dumz,zetap,dzeta,dzetab
      complex*16 f(10),ce(11),zr(10)
      integer i,j,nsol,nz,loop
      integer soln(10)
      logical polish
      logical polish_only,first_n_minus_2_roots_order_changed
      tol = 1.d-4

c The following are the original definition of coordinate
CC    setting the lens positions z1,z2,z3 and the source position zeta
CC    on the complex plane
      x1=da/2.
      x3=x1-db*dcos(phi/RADIAN)
      y3=db*dsin(phi/RADIAN)
      z1=dcmplx(-x1,0.)
      z2=dcmplx(x1,0.)
      z3=dcmplx(x3,y3)
      zeta=dcmplx(xs,ys)
      zetab=dcmplx(xs,-ys)
CC    converting the relative mass ratio to the total mass fraction
      q1=1./(1.+ma+mb)
      q2=ma/(1.+ma+mb)
      q3=mb/(1.+ma+mb)

!! Define the origin of coordinate system on the central star
c        x3 = -db*dcos(phi/RADIAN)
c        y3 = db*dsin(phi/RADIAN)
c        z1 = dcmplx(0.d0,0.d0)
c        z2 = dcmplx(-da,0.d0)
c        z3 = dcmplx(x3,y3)
c      zeta = dcmplx(xs,ys)
c      zetab = dcmplx(xs,-ys)
c! the mass must be normalised ?
c        q1 = 1.d0/(1.+ma+mb)
c        q2 = ma/(1.+ma+mb)
c        q3 = mb/(1.+ma+mb)

CC    some intermediate quantities using for the polynomial coefficient
      a=z1      +z2      +z3
      b=   z2*z3   +z3*z1   +z1*z2
      c=z1*z2*z3
      r=z1   *q1+z2   *q2+z3   *q3
      s=z1*z1*q1+z2*z2*q2+z3*z3*q3
      p=a-r
      q=b-r*a+s
      call conj3(ab,a)
      call conj3(bb,b)
      call conj3(cb,c)
      call conj3(pb,p)
      call conj3(qb,q)

      alpha=ab   -3*zetab
      beta =bb-2*ab*zetab +3*zetab**2
      gamma=cb  -bb*zetab+ab*zetab**2-zetab**3
      delta=pb   -2*zetab
      eta  =qb  -pb*zetab   +zetab**2

CC    f1=(z-z1)(z-z2)(z-z3) and f2=q1(z-z2)(z-z3)+q2(z-z3)(z-z1)+q3(z-z1)(z-z2)
CC    polynomial coefficients of f1^3
      f(1)=-c*c*c
      f(2)= 3*b*c*c
      f(3)=-3*c*(a*c+b*b)
      f(4)= 6*a*b*c+b*b*b+3*c*c
      f(5)=-3*(a*b*b+2*b*c+a*a*c)
      f(6)= 3*(2*a*c+a*a*b+b*b)
      f(7)=-(6*a*b+3*c+a*a*a)
      f(8)= 3*(a*a+b)
      f(9)=-3*a
      f(10)=dcmplx(1.,0.)

CC    the first part of the coefficient of the polynomialized lens equation
      ce(1) =              (eta-zeta*gamma)*f(1)
      do 10 i=2,10
      j = i - 1
      ce(i) = gamma*f(j) + (eta-zeta*gamma)*f(i)
 10   continue
      ce(11)= gamma*f(10)

CC    polynomial coefficients of f1^2*f2
      f(1)= c*q*c
      f(2)=-(2*b*c*q+c*c*p)
      f(3)= 2*c*(a*q+b*p)+b*b*q+c*c
      f(4)=-(2*q*(a*b+c)+2*c*(b+a*p)+b*b*p)
      f(5)= 2*p*(a*b+c)+2*(a*c+b*q)+a*a*q+b*b
      f(6)=-(2*(a*b+c)+2*(b*p+a*q)+a*a*p)
      f(7)= 2*a*p+2*b+q+a*a
      f(8)=-(a+p+a)
      f(9)= dcmplx(1.,0.)

CC    the second part of the coefficient of the polynomialized lens equation
      ce(1) = ce(1)             + (zeta*beta-delta)*f(1)
      do 20 i=2,9
      j = i - 1
      ce(i) = ce(i) - beta*f(j) + (zeta*beta-delta)*f(i)
 20   continue
      ce(10)= ce(10)- beta*f(9)

CC    polynomial coefficients of f1*f2^2
      f(1)=-c*q*q
      f(2)= 2*c*p*q+b*q*q
      f(3)=-(2*q*(b*p+c)+c*p*p+a*q*q)
      f(4)= 2*(b*q+c*p+a*p*q)+b*p*p+q*q
      f(5)=-(2*(a*q+b*p+p*q)+c+a*p*p)
      f(6)= 2*(a*p+q)+b+p*p
      f(7)=-(a+p+p)
      f(8)= dcmplx(1.,0.)

CC    the third part of the coefficient of the polynomialized lens equation
      ce(1) = ce(1)              + (1.-zeta*alpha)*f(1)
      do 30 i=2,8
      j = i - 1
      ce(i) = ce(i) + alpha*f(j) + (1.-zeta*alpha)*f(i)
 30   continue
      ce(9) = ce(9) + alpha*f(8)

CC    polynomial coefficients of f2^3
      f(1)= q*q*q
      f(2)=-3*p*q*q
      f(3)= 3*(p*p+q)*q
      f(4)=-(p*p*p+6*p*q)
      f(5)= 3*(p*p+q)
      f(6)=-3*p
      f(7)= dcmplx(1.,0.)

CC    the fourth part of the coefficient of the polynomialized lens equation
      ce(1) = ce(1)        + zeta*f(1)
      do 40 i=2,7
      j = i - 1
      ce(i) = ce(i) - f(j) + zeta*f(i)
 40   continue
      ce(8) = ce(8) - f(7)

CC    to solve the polynomialized lens equation using laguer3re method
CC    (see Press et al. - Numerical Recipe)
      nz = 10
      polish = .true.
      polish_only = .false.
c      call zroots3(ce,nz,zr,polish)
      call cmplx_roots_n
     * (zr, first_n_minus_2_roots_order_changed, ce, polish_only)
c     * (roots, first_n_minus_2_roots_order_changed, poly, polish_only)

CC    to check whether the solution of the polynomial is
CC    a solution of the (rational-form) lens equation 
 49   continue
      nsol = 0
      do 50 loop=1,10
         dumz = zr(loop)
         call lenseqt(zetap,dumz,z1,z2,z3,q1,q2,q3)
         err = zeta - zetap
         call dmag(derr,err)
         if(derr.gt.tol) then
            soln(loop) = 0
         else
            nsol = nsol + 1
            soln(loop) = 1
         endif
 50   continue

CC    to count the number of the solution, i.e. image and to check
CC    whether it is 4,6,8,10 (N+1 to N^2+1 with the interval of 2)
CC    otherwise to give error message
      if((nsol.lt.4).or.((nsol/2-int(nsol/2)).ne.0))then
ccc         write(6,*)'Wrong nsol', nsol
         do 51 loop=1,10
            dumz = zr(loop)
            call lenseqt(zetap,dumz,z1,z2,z3,q1,q2,q3)
            err = zeta - zetap
            call dmag(derr,err)
ccc            write(6,*)loop,derr
 51      continue
ccc         write(6,*)tol
c         write(6,*)'enter new tol'
c         read(5,*)tol
         tol = tol*9.5d-1
         if(tol.ge.4.d-3)go to 49
      endif

CC    to find the magnification for each image
CC    by calculating the inverse of the Jacobian at the image
      do 60 loop=1,10
         if(soln(loop).eq.1)then
            dumz = zr(loop)
            xi(loop) = dreal(dumz)
            yi(loop) = dimag(dumz)
            call lenseqtJ(dzeta,dumz,z1,z2,z3,q1,q2,q3)
            call conj3(dzetab,dzeta)
            detJ = 1. - dreal(dzeta*dzetab)
            am(loop) = 1./detJ
         else
            xi(loop) = 0.
            yi(loop) = 0.
            am(loop) = 0.
         endif
 60   continue

CC    to find the total magnification by adding the individual magnification.
      magtot = 0.
      do 61 loop=1,10
         magtot = magtot + dabs(am(loop))
 61   continue

      return
      end
CCCCC
      subroutine dmag(dm,z)
      complex*16 z
      double precision dm
      dm = dabs(dreal(z)) + dabs(dimag(z))
      return
      end
CCCCC
      subroutine conj3(zb,z)
      complex*16 zb,z
      double precision a,b
      a = dreal(z)
      b = dimag(z)
      zb = dcmplx(a,-b)
      return
      end
CCCCC
      subroutine lenseqt(lnt,z,z1,z2,z3,m1,m2,m3)
      complex*16 lnt,z,z1,z2,z3
      complex*16 zb,z1b,z2b,z3b
      double precision m1,m2,m3
      call conj3(z1b,z1)
      call conj3(z2b,z2)
      call conj3(z3b,z3)
      call conj3(zb,z)
      lnt = z + m1/(z1b-zb) + m2/(z2b-zb) + m3/(z3b-zb)
      return
      end
CCCCC
      subroutine lenseqtJ(lntJ,z,z1,z2,z3,m1,m2,m3)
      complex*16 lntJ,z,z1,z2,z3
      complex*16 zb,z1b,z2b,z3b
      double precision m1,m2,m3
      call conj3(z1b,z1)
      call conj3(z2b,z2)
      call conj3(z3b,z3)
      call conj3(zb,z)
      lntJ = m1/(z1b-zb)**2 + m2/(z2b-zb)**2 + m3/(z3b-zb)**2
      return
      end
CCCCC
      subroutine zroots3(a,m,roots,polish)
CU    uses laguer3
      integer MAXM
      real*8 EPS
      parameter (MAXM=15)
      parameter (EPS=1.d-14)
      integer m
      complex*16 a(m+1),roots(m)
      complex*16 ad(MAXM)
      complex*16 x,b,c
      integer i,j,jj,its
      logical polish

      do 10 j=1,m+1
         ad(j) = a(j)
 10   continue

        zero=0.d0

      do 12 j=m,1,-1
         x = cmplx(0.,0.)
         call laguer3(ad,j,x,its)
         if(abs(dimag(x)).le.2.*EPS**2*abs(dreal(x)))then
            x = cmplx(dreal(x),zero)
         endif
         roots(j) = x
         b = ad(j+1)
         do 11 jj=j,1,-1
            c = ad(jj)
            ad(jj) = b
            b = x*b + c
 11      continue
 12   continue

      do 13 j=1,m+1
         ad(j) = a(j)
 13   continue

      if(polish)then
         do 14 j=1,m
            call laguer3(ad,m,roots(j),its)
 14      continue
      endif

      do 17 j=2,m
         x = roots(j)
         do 15 i=j-1,1,-1
            if(dreal(roots(i)).le.dreal(x))go to 16
            roots(i+1) = roots(i)
 15      continue
         i = 0
 16      roots(i+1) = x
 17   continue

      return

      end


      subroutine laguer3(ad,m,x,its)
      integer MAXM,MAXIT,MR,MT
      real*8 EPSS
c      parameter (MAXM=15,MR=8,MT=10000,MAXIT=MT*MR)
      parameter (MAXM=15,MR=8,MT=100,MAXIT=MT*MR)
      parameter (EPSS=1.d-14)
c      parameter (EPSS=1.d-5)
      integer m,its
      complex*16 a(m+1),x
      integer i,iter,j
      real*8 abx,abp,abm,err,frac(MR),az,bz
      complex*16 ad(MAXM),dx,x1,b,d,f,g,h,gp,gm,cz
      save frac
      data frac /.5d0,.25d0,.75d0,.13d0,.38d0,.62d0,.88d0,1.d0/

      do 10 i=1,m+1
         a(i) = ad(i)
 10   continue

      do 12 iter=1,MAXIT
         its = iter
         b = a(m+1)
         err = abs(b)
         d = cmplx(0.,0.)
         f = cmplx(0.,0.)
         abx = abs(x)
         do 11 j=m,1,-1
            f = x*f + d
            d = x*d + b
            b = x*b + a(j)
            err = abs(b) + abx*err
 11      continue
         err=EPSS*err 
         if(abs(b).le.err)then
            return
         else
            g = d/b
            h = g*g - f/b - f/b
            gp = g + sqrt(m*m*h+g*g-m*h-m*g*g)
            gm = g - sqrt(m*m*h+g*g-m*h-m*g*g)
            abp = abs(gp)
            abm = abs(gm)
            if(abp.lt.abm)gp=gm
            if(max(abp,abm).gt.0.)then
               dx = m/gp
            else
               az=dlog(1.+abx)
               bz=real(iter)
               cz=cmplx(az,bz)
               dx = exp(cz)
            endif
         endif
         x1 = x - dx
         if(x.eq.x1)return
         if(mod(iter,MT).ne.0)then
            x = x1
         else
            x = x - dx*frac(iter/MT)
         endif
 12   continue

c      write(6,*) 'too many iterations in SUBROUTINE laguer3'

      return

      end

