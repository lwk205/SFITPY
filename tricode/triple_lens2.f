      subroutine GETMAG
     +(da, db, phi, ma, mb,
     + scx, scy, rho, magnif, lgamma,
     + mgridsize, xorigin, yorigin,
     + errorflag,zr)
      implicit none
      real*8 sx, sy, magnif, sx_temp, sy_temp, rho2_temp, s0
      integer ii
      integer SPMAX, SP, dimen
      parameter(SPMAX = 500)
      real*8 POLYFACTOR
      parameter(POLYFACTOR = 1.5d0)
      integer max_sol
      parameter(max_sol = 10)
      integer mparity(max_sol)
      real*8 rhopoly, rho2, rho, lgamma, scx, scy, drho2, rhopoly2
      real*8 mgridsize, xorigin, yorigin
      real*8 q1, q2, q3, ma, mb, da, db, phi, x1, x2, x3, y3
      real*8 maxdis
      integer i, j, k, z
      integer errorflag, ctrl
      logical verbose
      real*8 PI
      parameter (PI=3.141592653589793238462643d0)
      real*8 RADIAN
      parameter(RADIAN=1.8d2/PI)
      real*8 ix(SPMAX+1, max_sol), iy(SPMAX+1,max_sol)
      integer*4 nsol(SPMAX+1),parity(SPMAX+1,max_sol)
      real*8 ix_temp(SPMAX+1, max_sol), iy_temp(SPMAX+1,max_sol)
      integer*4 nsol_temp(SPMAX+1),parity_temp(SPMAX+1,max_sol)
      integer i1, max_n, max_i, min_n, msol
      integer ileft, jleft, iright, jright
      integer*4 CNUM
      integer*4 BDTEST
      parameter (BDTEST=100)
      parameter (CNUM=100000)
      real*8 magtott(max_sol),amt(max_sol),ixt(max_sol),iyt(max_sol)
      real*8 swapx(max_sol),swapy(max_sol),swapp(max_sol)
      integer m,n
      integer*4 head(SPMAX),tail(SPMAX),headtail(SPMAX)
      real*8 loopx(SPMAX*max_sol+max_sol,max_sol),
     + loopy(SPMAX*max_sol+max_sol,max_sol)
      real*8 vectorx(SPMAX*max_sol+max_sol,max_sol),
     + vectory(SPMAX*max_sol+max_sol,max_sol),sarea
      integer*4 loopp(max_sol),loopend(max_sol)
      integer*4 minim
      integer*4 tag, cstart, lstart, istart,startl,endl
      integer*4 indexs,indexe,initl,stringc
      integer*4 gridnumx,gridnumy
      real*8 gridminx,gridminy,gridx,gridy
      real*8 miniix,maxix,miniiy,maxiy
      real*8 loopminiix(max_sol),loopmaxix(max_sol)
      real*8 loopminiiy(max_sol),loopmaxiy(max_sol)
      real*8 sxi,syi
      real*8 marea
      real*8 gin(max_sol*SPMAX),gout(max_sol*SPMAX)
      integer*4 m3,m4,m5,m6,m7,mmax,mmin,mint,mindex,m8,mmmin
      integer*4 segp(100*SPMAX+100)
      real*8 dis,seglim,mmid
      real*8 segmax(100*SPMAX+100)
      real*8 segmin(100*SPMAX+200)
      real*8 segx1(100*SPMAX+100)
      real*8 segk(100*SPMAX+100)
      real*8 segkk,deltax,segxx1,segmmin,segmmax
      real*8 seglen2(max_sol*SPMAX+max_sol,max_sol)
      integer*4 mpoint(100*SPMAX+100),mindex2
      integer*4 segpp
      integer*4 check
      integer*4 gstep
      real*8 mz
      integer*4 sindex(max_sol)
      integer*4  nsolpos,nsolneg,insol
      logical bad(SPMAX+1)
      integer SP_count, k0, k1, k2, k3, k4, k5 
      integer ncross, errortype
      parameter (ncross = 5)
      real*8 accuracy, dis_max

      integer*4 ntrys
      complex*16 zr(10)

      errorflag=0
      verbose  =.false.
      ctrl     =1
      ntrys =0 !number of times through the 15 continue loop

      rhopoly = rho*POLYFACTOR
      rho2  = rho*rho

      q1 = 1.d0/(1.+ma+mb)
      q2=ma/(1.+ma+mb)
      q3=mb/(1.+ma+mb)
cc Below is the initial geometry
      x1=-da/2.d0
      x2=da/2.d0
      x3=x1-db*dcos(phi/RADIAN)
      y3=db*dsin(phi/RADIAN)
c Now begins my geometry
c      x1 = 0.d0
c      x2 = -da
c      x3 = -db*dcos(phi/RADIAN)
c      y3 = db*dsin(phi/RADIAN)
c     Define Necessary Quantities for Triple Lens Ray-Shooting

      SP    = SPMAX
c     NOTE: SP !(necessarily)= SPMAX
      dimen = SPMAX

      maxdis = 1d0/SP*4.d0 
c     Related to RayShooting
c     change

15    continue
      k = 1
      z = 1
c     initializing k & z

      ntrys=ntrys+1
c     Avoid infinite loop scenerio
      if (ntrys .gt. 100) then
         errorflag=1
c         write(6,*) ntrys,'>=100'
         go to 1000
      endif

      do 500 j=1,max_sol

      loopp(j)   = 0.d0
      loopend(j) = 0.d0      

      do 510 i=1,SPMAX+1
      parity(i,j) = 0.d0
      parity_temp(i,j) = 0.d0
      bad(i) = .true.
510   continue

500   continue

      i1 = 0

      max_n = 0
      max_i = 0
      min_n = max_sol + 1

      accuracy = rhopoly*dsin(PI/SPMAX)*2.d0

      do 10 i=1,SPMAX

      sx=scx+rhopoly*dcos(2d0*PI/SPMAX*(i-1))
      sy=scy+rhopoly*dsin(2d0*PI/SPMAX*(i-1))

c      call getimt(ixt,iyt,amt,msol,da,db,phi,ma,mb,sx,sy,magtott)
      call lenssolver(ixt,iyt,amt,msol,da,db,phi,ma,mb,sx,sy,magtott,
     + mparity,max_sol,errortype, dis_max,
     + x1, x2, x3, y3, q1, q2, q3,zr)
      if(errortype.ne.0.or.dis_max.gt.accuracy) then
      do j = 1,ncross - 1
      sx=scx+rhopoly*dcos(2d0*PI/SPMAX*(i-1) + 2d0*PI/SPMAX/ncross*j)
      sy=scy+rhopoly*dsin(2d0*PI/SPMAX*(i-1) + 2d0*PI/SPMAX/ncross*j) 
      call lenssolver(ixt,iyt,amt,msol,da,db,phi,ma,mb,sx,sy,magtott,
     + mparity,max_sol,errortype, dis_max,
     + x1, x2, x3, y3, q1, q2, q3,zr)
      if(errortype.eq.0.and.dis_max.le.accuracy) then
      goto 713
      endif
      enddo
      rhopoly = rhopoly*POLYFACTOR
      goto 15
      endif
  
713   continue

      i1 = i1 + 1
      nsol_temp(i1) = msol

      if(msol.gt.max_n) then
      max_n = msol
      max_i = i1
      endif

      if(msol.lt.min_n) then
      min_n = msol
      endif

      do 31 j=1,max_sol
      ix_temp(i1,j)    = ixt(j)
      iy_temp(i1,j)    = iyt(j)
      parity_temp(i1,j)= mparity(j)
31    continue

11    continue
10    continue

      SP = i1 

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(max_n.eq.min_n) then

      do i = 1,SP
      nsol(i)     = nsol_temp(i)
      do j = 1,max_sol
      ix(i,j)     = ix_temp(i,j)
      iy(i,j)     = iy_temp(i,j)
      parity(i,j) = parity_temp(i,j)
      enddo
      enddo
      goto 12
      endif

      ileft  = max_i
      iright = max_i

      do 13 j = 1,SP-1
      ileft = ileft - 1
      if(ileft.lt.1) then
      jleft = ileft + SP
      else
      jleft = ileft
      endif
      if(nsol_temp(jleft).lt.max_n) goto 21
13    continue

21    continue

      do 14 j = 1,SP-1
      iright = iright + 1
      if(iright.gt.SP) then
      jright = iright - SP
      else
      jright = iright
      endif
      if(nsol_temp(jright).lt.max_n) goto 22
14    continue
22    continue

      s0 = (ileft + iright)/2
      if(s0.lt.1)  s0 = s0 + SP
      if(s0.gt.SP) s0 = s0 - SP

      
      do i = 1,SP
      ii = s0 + i - 1
      if(ii.gt.SP) ii = ii - SP
      nsol(i)     = nsol_temp(ii)
      do j = 1,max_sol
      ix(i,j)     = ix_temp(ii,j)
      iy(i,j)     = iy_temp(ii,j)
      parity(i,j) = parity_temp(ii,j)
      enddo
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c To Move The Element With Maximum Number of Sol. in the Beginning 
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
12    continue

      goto 888
      do i = 1,SP
      nsol_temp(i) = nsol(i)
      bad(i) = .false.
      do j=1,max_sol
      ix_temp(i,j) = ix(i,j)
      iy_temp(i,j) = iy(i,j)
      parity_temp(i,j) = parity(i,j)
      enddo
      enddo

      do i = 1,SP
      j = i + 1
      if(j.gt.SP) j = j - SP
      if(nsol(i).ne.nsol(j)) then
c      k0 = i - 2
c      k1 = i - 1
      k2 = i
      k3 = i + 1
c      k4 = i + 2
c      k5 = i + 3
c      if(k0.lt.1)  k0 = k0 + SP
c      if(k1.lt.1)  k1 = k1 + SP
      if(k3.gt.SP) k3 = k3 - SP
c      if(k4.gt.SP) k4 = k4 - SP
c      if(k5.gt.SP) k5 = k5 - SP
c      bad(k0) = .true.
c      bad(k1) = .true.
      bad(k2) = .true.
      bad(k3) = .true.
c      bad(k4) = .true.
c      bad(k5) = .true.
      endif
      enddo

      SP_count = 0

      do 111 i = 1,SP
      if(bad(i)) then
      goto 111
      else
      SP_count = SP_count + 1
      nsol(SP_count) = nsol_temp(i)
      do j=1,max_sol
      ix(SP_count,j) = ix_temp(i,j)
      iy(SP_count,j) = iy_temp(i,j)
      parity(SP_count,j) = parity_temp(i,j)
      enddo
      endif
111   continue

      SP = SP_count

888   continue

      nsol(SP + 1) = nsol(1)
      do j = 1,max_sol
      ix(SP + 1,j) = ix(1,j)
      iy(SP + 1,j) = iy(1,j)
      parity(SP + 1,j) = parity(1,j)
      enddo

cc      do j=1,max_sol
cc      do i=1,SP+1
cc      write(70+j,801) ix(i,j), iy(i,j), parity(i,j)
cc      enddo
cc      close(70+j)
cc      enddo
cc701   format(2f12.7, i5)



ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do 45 j=1,max_sol
      if (parity(1,j).NE.0) then
      head(k)=1
      head(k+1)=j
      k=k+2
      endif
45    continue

      do 40 i=1,SP

      do 60 j=1,max_sol
      swapx(j)=0
      swapy(j)=0
      swapp(j)=0
60    continue

      if (nsol(i+1).EQ.nsol(i)) then
      call COMP2(ix,iy,parity,i,i+1,dimen,sindex,max_sol)

      do 50 j=1,max_sol
      swapx(j)=ix(i+1,sindex(j))
      swapy(j)=iy(i+1,sindex(j))
      swapp(j)=parity(i+1,sindex(j))
50    continue

      do 51 j=1,max_sol
      ix(i+1,j)=swapx(j)
      iy(i+1,j)=swapy(j)
      parity(i+1,j)=swapp(j)
51    continue

      else if (nsol(i+1).GT.nsol(i)) then

      call COMP2(ix,iy,parity,i,i+1,dimen,sindex,max_sol)



      do 81 j=1,max_sol
      swapx(j)=ix(i+1,sindex(j))
      swapy(j)=iy(i+1,sindex(j))
      swapp(j)=parity(i+1,sindex(j))
81    continue
 
      do 82 j=1,max_sol
      ix(i+1,j)=swapx(j)
      iy(i+1,j)=swapy(j)
      parity(i+1,j)=swapp(j)
82    continue

      do 80 j=1,max_sol

      if (parity(i,j).EQ.0.AND.parity(i+1,j).NE.0) then
           head(k)=i+1
           head(k+1)=j
           k=k+2
      endif
80    continue

      else

      call COMP2(ix,iy,parity,i+1,i,dimen,sindex,max_sol)



         do 90 j=1,max_sol
         swapx(sindex(j))=ix(i+1,j)
         swapy(sindex(j))=iy(i+1,j)
         swapp(sindex(j))=parity(i+1,j)
90       continue

         do 91 j=1,max_sol
         ix(i+1,j)=swapx(j)
         iy(i+1,j)=swapy(j)
         parity(i+1,j)=swapp(j)
91       continue

         do 100 j=1,max_sol
         if (swapp(j).EQ.0.AND.parity(i,j).NE.0) then
         tail(z)=i
         tail(z+1)=j
         z=z+2
         endif
100      continue

      endif

40    continue

      do 110 j=1,max_sol

      if ((parity(SP+1,j).NE.0)) then
      tail(z)=SP+1
      tail(z+1)=j
      z=z+2
      endif


110   continue

      k=k-1
      z=z-1

ccccccccccccc link the loop cccccccccccccccccccccccccccccccccccccccc

      n=1
      do 170 i=1,(k/2)
         minim=SP+1
         do 180 j=1,(z/2)

         if (head(i*2).EQ.tail(j*2)) then
            if (tail(j*2-1).GE.head(i*2-1)) then
               if ((tail(j*2-1)-head(i*2-1)).LT.minim) then
         minim=tail(j*2)-head(i*2)
         headtail(n)=i
         headtail(n+1)=j
               endif
            endif
        endif

180   continue
      n=n+2
170   continue

      n=n-1

160    continue

      do 161 j=1,max_sol
         loopp(j)=0
         do 162 i=1,SP*max_sol+max_sol
         loopx(i,j)=0d0
         loopy(i,j)=0d0
162   continue
161   continue

      tag=0

      do 190 k=1,max_sol
          do 200 z=1,(n/2)

             if(headtail(2*z-1).NE.0) then
             cstart=head(2*headtail(2*z-1))
             lstart=head(2*headtail(2*z-1)-1)
             istart=2*z-1
             loopp(k)=parity(lstart,cstart)
             goto 210
             endif

200        continue

        initl=1
        goto 195

210     indexs=istart
        indexe=istart+1
        initl=1

        do 220 z=1,n
        
        if((indexe-indexe/2*2).EQ.0) then

        tag=-1

        startl=head(2*headtail(indexs)-1)
        stringc=head(2*headtail(indexs))
        endl=tail(2*headtail(indexe)-1)
                        
        else

        tag=1

        startl=tail(2*headtail(indexs)-1)
        stringc=tail(2*headtail(indexs))
        endl=head(2*headtail(indexe)-1)

        endif

      call cpstring(ix,iy,loopx,loopy,dimen,
     +                startl,endl,stringc,initl,k,max_sol)

       headtail(indexs)=0
       headtail(indexe)=0

      call findp (headtail,head,tail,ix,iy,
     +            parity,dimen,n,loopp(k),
     +            stringc,endl,cstart,lstart,indexs,indexe,max_sol)

      if(indexs.EQ.0) then
      goto 195

      endif

220    continue

195    loopend(k)=initl-1

190    continue

      do 230 i=1,n

      if(headtail(i).NE.0) then
      ctrl=ctrl+1
      rhopoly=rhopoly*POLYFACTOR
      if(verbose) then
      write(6,*) 'Error 2.1 Head-Tail Fault'
      write(6,*) i,n,rhopoly
      endif

    
      if(ctrl.LT.CNUM) then 
         goto 15
         else 
         write(6,*) 'Error 2.0 Head Tail Errors Reach Critical Level'
         write(6,*) 'Point Calculation Returned'
         errorflag=2
         goto 1000
         endif

      endif     

230   continue


      
      z=0

      do 235 j=1,max_sol

cc      write(6,*) j, loopend(j)
     
      if (loopend(j).NE.0) then
      
      z=z+1
  
         if (loopp(j).LT.0)   then

      call reversearray(loopx,j,loopend(j),(dimen*max_sol+max_sol)
     + ,max_sol)
      call reversearray(loopy,j,loopend(j),(dimen*max_sol+max_sol)
     + ,max_sol)
         endif

cc      write(6,*) 'write loop'

c      do i=1,loopend(j)
c      write(80+j,801) loopx(i,j), loopy(i,j)
c      enddo
c801   format(2f12.7, i5)
c      close(80+j)

      endif

235   continue

c      write(6,*) rhopoly, mgridsize
c      pause
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c The inverse shoot                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Boundary Check                                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      mgridsize=0.00004d0*(fstol/0.001d0)**(2d0/3d0)*(rho/0.001d0)
cc      mgridsize=(fstol*sqrt(6d0*PI*magp))**(2d0/3d0)*rho

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      gstep=BDTEST
      
      do 400 j=1,max_sol

      if (loopend(j).NE.0) then
   
      do 410 i=1,loopend(j)-1
      call boundarytest(check,loopx(i,j),loopy(i,j)
     +               ,loopx(i+1,j),loopy(i+1,j),gstep,
     +               rho2,scx,scy,
     + x1, x2, x3, y3, q1, q2, q3)
      
      if (check.EQ.0) then
      rhopoly=rhopoly*POLYFACTOR
      ctrl=ctrl+1
      if(verbose) then
      write(6,*) 'Error 3.2 Boundary Test Fails'
      write(6,*) rhopoly,SP
      endif
      if(ctrl.LT.CNUM) then
         goto 15
         else 
         write(6,*) 'Error 3.0 Boundary Errors Reach Critical Level'
         errorflag=3
         goto 1000
         endif

      endif

410   continue
      endif
400   continue 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 600 j=1,max_sol
      if(loopend(j).NE.0) then
      loopminiiy(j)=loopy(1,j)
      loopmaxiy(j)=loopy(1,j)
      loopminiix(j)=loopx(1,j)
      loopmaxix(j)=loopx(1,j)

      do 610 i=1,loopend(j)-1
      vectorx(i,j)=loopx(i+1,j)-loopx(i,j)
      vectory(i,j)=loopy(i+1,j)-loopy(i,j)
      seglen2(i,j)=vectorx(i,j)*vectorx(i,j)+
     +             vectory(i,j)*vectory(i,j)
      if(loopy(i,j).GT.loopmaxiy(j)) then
      loopmaxiy(j)=loopy(i,j)
      endif
      if(loopx(i,j).GT.loopmaxix(j)) then
      loopmaxix(j)=loopx(i,j)
      endif
      if(loopy(i,j).LT.loopminiiy(j)) then
      loopminiiy(j)=loopy(i,j)
      endif
      if(loopx(i,j).LT.loopminiix(j)) then
      loopminiix(j)=loopx(i,j)
      endif
610   continue
      vectorx(loopend(j),j)=vectorx(1,j)
      vectory(loopend(j),j)=vectory(1,j)
      seglen2(loopend(j),j)=seglen2(1,j)

      if((loopmaxix(j)-loopminiix(j)).LT.mgridsize) then
         if((loopmaxiy(j)-loopminiiy(j)).LT.mgridsize) then
      call stokes(loopx,loopy,vectorx,vectory,loopend(j),j,sarea,
     +            dimen,max_sol)
           
cccccc           write(6,*) 'small area detected',j
           if(sarea.LE.0) then
           loopend(j)=0
           endif
         endif
      endif
      
      endif
600   continue
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mz=0d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c      mgridsize=0.00004d0*(fstol/0.001d0)**(2d0/3d0)*(rho/0.001d0)
cc     mgridsize=(fstol*sqrt(6d0*PI*magp))**(2d0/3d0)*rho
cc      yorigin=0d0-0.5d0*mgridsize
cc      xorigin=0d0-offset
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      m7=0
      do 921 j=1,max_sol
         if(loopend(j).NE.0) then
          maxiy=loopy(1,j)
          miniiy=loopy(1,j)
         endif
921   continue
      do 900 j=1,max_sol
         if(loopend(j).NE.0) then
      do 910 i=1,loopend(j)-1
      
      if(maxiy.LT.loopy(i,j)) then
      maxiy=loopy(i,j)
      endif
      if(miniiy.GT.loopy(i,j)) then
      miniiy=loopy(i,j)
      endif
      dis=loopy(i+1,j)-loopy(i,j) 
      segkk=(loopx(i+1,j)-loopx(i,j))/(dis)
       
      if(dis.GT.0) then
      segmmin=loopy(i,j)
      segmmax=loopy(i+1,j)
      segxx1=loopx(i,j)
      segpp=1
      else if (dis.LT.0) then
      segmmin=loopy(i+1,j)
      segmmax=loopy(i,j)
      segxx1=loopx(i+1,j)
      dis=-dis
      segpp=-1
      else
      goto 910
      endif
      
      if(dis.LE.maxdis) then
      
      m7=m7+1
      segmax(m7)=segmmax
      segmin(m7)=segmmin
      segp(m7)=segpp
      segx1(m7)=segxx1
      segk(m7)=segkk
      
      else
      
      deltax=segkk*maxdis
      
      do 951 m8=0,nint(dis/maxdis-0.5d0)-1
      m7=m7+1
      segmin(m7)=segmmin+m8*maxdis
      segmax(m7)=segmin(m7)+maxdis
      segp(m7)=segpp
      segk(m7)=segkk
      segx1(m7)=segxx1+deltax*m8
951   continue
      m7=m7+1
      segmin(m7)=segmax(m7-1)
      segmax(m7)=segmmax
      segp(m7)=segpp
      segk(m7)=segkk
      segx1(m7)=segx1(m7-1)+deltax
      
      endif
910   continue
         endif
900   continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 919 m8=1,m7
      mpoint(m8)=m8
919   continue
      
c      write(6,*) 'm7=',m7

      call hpsorti(m7,segmin,mpoint)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c      DELAPSE =DTIME(TIMEARRAY)
c      write(6,*) 'sort',DELAPSE,'Poly:',SP,'maxdis:',maxdis
c      DELAPSE =DTIME(TIMEARRAY)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call buildgrid(miniiy,maxiy,yorigin,mgridsize,
     +          gridminy,gridnumy)
      
      gridy=gridminy
      mmmin=0

      do 630 i=1,gridnumy

      m3=0
      m4=0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   binary search                                               c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call  locate(segmin,mmmin,m7,gridy,mmin)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mmmin=mmin-1
      mindex=mmin
      mindex2=mpoint(mindex)
      seglim=segmin(mmin)-maxdis
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c!!!!!If mindex = 1, this gives an array overrun, so I changed the condition
c!!!!!to mindex>=2. ad hoc. JCY
c990   if(mindex.GE.1) then
990   if(mindex.GE.2) then
       if(seglim.LE.segmin(mindex)) then
         if(gridy.LT.segmax(mindex2)) then
           if(segp(mindex2).EQ.1) then
              m4=m4+1
            gout(m4)=segk(mindex2)*(gridy-segmin(mindex))+segx1(mindex2)
           else
              m3=m3+1
            gin(m3)=segk(mindex2)*(gridy-segmin(mindex))+segx1(mindex2)
           endif
        endif
        mindex=mindex-1
        mindex2=mpoint(mindex)
        goto 990 
      endif
      endif

      if(m3.NE.0) then   
      
      if(m3.NE.m4) then
      if(verbose) then
      write(6,*) 'Error 4.1 Number of Enter & Exit Points Do Not Match'
      write(6,*) 'm3=',m3,'m4=',m4,'rhopoly=',rhopoly,'gridy=',gridy
      endif
      rhopoly=rhopoly*POLYFACTOR
      ctrl=ctrl+1
      if(ctrl.LT.CNUM) then 
         goto 15
         else 
         write(6,*) 'Error 4.0 Enter & Exit Number Match Errors Reach 
     +               Critical Level'
         errorflag=4
         goto 1000
         endif
      endif

      call sort(gin,m3)
      call sort(gout,m4)
      
      do 650 m5=1,m3
      if(gin(m5).GT.gout(m5)) then
      if(verbose) then
      write(6,*) 'Error 5.2 Enter > Exit (Possible Crossings)'
      write(6,*) 'rhopoly=',rhopoly,gin(m5),gout(m5),gridy
      endif
      rhopoly=rhopoly*POLYFACTOR
      ctrl=ctrl+1
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      open(13, FILE='loops.dat', STATUS='unknown')
c13      write(13,*) z
c13      write(13,*) '-------------------------------------------'
c13      do 1260 jdebug=1,5
c13      if (loopend(jdebug).NE.0) then
c13      write(13,*) loopend(jdebug),loopp(jdebug)
c13      write(13,*) '-------------------------------------------'
c13      do 1270 idebug=1,loopend(jdebug)
c13      write(13,*) loopx(idebug,jdebug),loopy(idebug,jdebug)
c131270   continue
c13      write(13,*) '--------------------------------------------'
c13      endif
c131260   continue
c      close(13)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      if(ctrl.LT.CNUM) then 
         goto 15
         else 
         write(6,*) 'Error 5.0 Enter > Exit Errors Reach Critical Level'
         errorflag=5
         goto 1000
         endif
      endif

      call buildgrid(gin(m5),gout(m5),xorigin,
     +               mgridsize,gridminx,gridnumx)

      gridx=gridminx

      do 660 m6=1,gridnumx

c      call inverseshoot(sxi,syi,gridx,gridy,bdiv2,u1,u2)
      call
     + invshoot(sxi, syi, gridx, gridy, x1, x2, x3, y3, q1, q2, q3)
c      write(92,902) gridx, gridy, sxi, syi

901   format(2f16.7)
902   format(4f16.7)


      drho2=(sxi-scx)**2+(syi-scy)**2
c      if(drho2.le.1e-6) then
c      write(93,902) gridx, gridy, sxi, syi
c      endif
      if (drho2.le.rho2) then
c       mz=mz+((1d0-lgamma)+3d0/2d0*lgamma*
c     &              dsqrt(1d0-(drho2/rho2))) 
c      write(91,902) gridx, gridy, sxi, syi
      if(lgamma.ne.0) then
      mz=mz+((1d0-lgamma)+3d0/2d0*lgamma*
     +   sqrt(1d0-drho2/rho2))
      else
      mz=mz+1d0
      endif
c      mz=mz+1
c      write(11,*) gridx,gridy
      endif

      gridx=gridx+mgridsize
660   continue
 
650   continue
 
      do 651 m5=1,m3-1
      if (gout(m5).GT.gin(m5+1)) then
      if(verbose) then
      write(6,*) 'Error 6.1 Exit > Next Enter'
      write(6,*) rhopoly,gout(m5),gin(m5+1),gridy
      endif
      rhopoly=rhopoly*POLYFACTOR
      ctrl=ctrl+1
      
       if(ctrl.LT.CNUM) then
         goto 15
         else
         write(6,*) 'Error 6.0 Exit > Next Enter Errors reach critical 
     +               value'
         errorflag=6
         goto 1000
       endif
      endif
651   continue   
      endif
      gridy=gridy+mgridsize
630   continue
c      DELAPSE =DTIME(TIMEARRAY)
c      write(6,*) 'grid determination (new method)',DELAPSE,'Poly:',SP
      marea=mz*mgridsize*mgridsize
      magnif=marea/(PI*rho2)

1000  return

      END

      subroutine findp (headtail,head,tail,ix,iy,
     +                  parity,dimen,n,loopp,
     +                  column,endl,cstart,lstart,indexs,indexe
     + ,max_sol)

      implicit none
      integer*4 max_sol
      integer*4 dimen,tag,column,endl,cstart,lstart,indexs,indexe
      integer*4 headtail(dimen),head(dimen),tail(dimen)
      integer*4 parity(dimen+1,max_sol),i,n,loopp,lpar,ll,lc
      real*8 ix(dimen+1,max_sol),iy(dimen+1,max_sol),mini
      real*8 x,y,ld
      
      x=ix(endl,column)
      y=iy(endl,column)
c      mini=(ix(lstart,cstart)-x)**2+(iy(lstart,cstart)-y)**2
      mini=10000
      indexs=0
      indexe=0
      if((ix(lstart,cstart)-x)**2+(iy(lstart,cstart)-y)**2.EQ.0) then
      if ((lstart.NE.endl).OR.(cstart.NE.column)) then
      goto 20
      endif
      endif

      do 10 i=1,n
         if(headtail(i).NE.0) then
           if((i-i/2*2).EQ.1) then
             ll=head(2*headtail(i)-1)
             lc=head(2*headtail(i))
             lpar=parity(ll,lc)
             if (lpar.EQ.(loopp)) then
                ld=(ix(ll,lc)-x)**2+(iy(ll,lc)-y)**2
                if (ld.LE.mini) then
                   mini=ld
                   indexs=i
                   indexe=i+1
                endif
             endif
           else
             ll=tail(2*headtail(i)-1)
             lc=tail(2*headtail(i))
             lpar=parity(ll,lc)
             if (lpar.EQ.((-1)*(loopp))) then
                ld=(ix(ll,lc)-x)**2+(iy(ll,lc)-y)**2
                if (ld.LE.mini) then
                   mini=ld
                   indexs=i
                   indexe=i-1
                endif
             endif
           endif
         endif
10    continue

20    return
      end

       subroutine cpstring(stringx1,stringy1,stringx2,stringy2,dimen,
     +                    startl1,endl1,c1,startl2,c2, max_sol)

      implicit none
      integer*4 max_sol
      integer*4 dimen,startl1,endl1,c1,startl2,c2,i
      real*8 stringx1(dimen+1,max_sol),stringy1(dimen+1,max_sol)
      real*8 stringx2((dimen+1)*max_sol,max_sol)
     + ,stringy2((dimen+1)*max_sol,max_sol)

      if(endl1.GT.startl1) then

      do 10 i=0,(endl1-startl1)

      stringx2(startl2+i,c2)=stringx1(startl1+i,c1)
      stringy2(startl2+i,c2)=stringy1(startl1+i,c1)

10    continue

      startl2=startl2+endl1-startl1+1

      else

      do 20 i=0,(startl1-endl1)

      stringx2(startl2+i,c2)=stringx1(startl1-i,c1)
      stringy2(startl2+i,c2)=stringy1(startl1-i,c1)

20    continue

      startl2=startl2+startl1-endl1+1

      endif

      return
      end

      subroutine COMP2(mx,my,mp,l1,l2,dimen,sindex, max_sol)

c     find the nearest points from line l1 to line l2, return an index matrix referring from the column number of l1 to that
c     l2

      implicit none
      integer*4 max_sol
      integer*4 dimen
      real*8 mx(dimen+1,max_sol),my(dimen+1,max_sol),mini,dist
      integer*4 mp(dimen+1,max_sol),l1,l2,i,j,k,sindex(max_sol),i2,j2
      integer*4 s1(max_sol),s2(max_sol)

      do k=1,max_sol
      sindex(k)=-1
      if(mp(l1,k).ne.0) then
      s1(k)=1
      else
      s1(k)=-1
      endif

      if(mp(l2,k).ne.0) then
      s2(k)=1
      else
      s2(k)=-1
      endif

      enddo

      i2=-1
      j2=-1

      do 30 k=1,max_sol
      mini=10000d0
      do 10 i=1,max_sol
        if(s1(i).eq.1) then
            do 20 j=1,max_sol
              if(s2(j).eq.1) then
                if(mp(l1,i).eq.mp(l2,j)) then
                dist=(mx(l1,i)-mx(l2,j))**2+(my(l1,i)-my(l2,j))**2
                  if(dist.LT.mini) then
                  mini=dist
                  i2=i
                  j2=j
                  sindex(i2)=j2
                  endif
                endif
              endif
20           continue
           endif
10     continue
       if((i2.ne.-1).and.(j2.ne.-1)) then
       s1(i2)=0
       s2(j2)=0
       endif
30     continue

      do 40 k=1,max_sol
      if(s1(k).eq.-1) then
         do 50 i=1,max_sol
         if(s2(i).ne.0) then
         sindex(k)=i
         s1(k)=0
         s2(i)=0
         goto 40
         endif
50    continue
      endif
40    continue

      return
      end


      
      subroutine COMP(mx,my,mp,l1,c1,l2,c2,dimen,max_sol)

      implicit none
      integer*4 max_sol
      integer*4 dimen
      real*8 mx(dimen+1,max_sol),my(dimen+1,max_sol),mini,dist
      integer*4 mp(dimen+1,max_sol),l1,c1,l2,c2,h

      mini=10000d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Consider Revise                                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do 10 h=1,max_sol

      if (mp(l2,h).EQ.mp(l1,c1)) then
         dist=(mx(l1,c1)-mx(l2,h))**2+(my(l1,c1)-my(l2,h))**2
         if (dist.LT.mini) then

      mini=dist
      c2=h

      endif
      endif

10    continue
      return
      end
       
      subroutine SWAP(mx,my,mp,l1,c1,l2,c2,dimen,max_sol)

      implicit none
      integer*4 max_sol
      integer*4 dimen
      real*8 mx(dimen+1,max_sol),my(dimen+1,max_sol)
      real*8 mswap
      integer*4 h,mp(dimen+1,max_sol),l1,c1,l2,c2

      mswap=mx(l1,c1)
      mx(l1,c1)=mx(l2,c2)
      mx(l2,c2)=mswap

      mswap=my(l1,c1)
      my(l1,c1)=my(l2,c2)
      my(l2,c2)=mswap

      mswap=mp(l1,c1)
      mp(l1,c1)=mp(l2,c2)
      mp(l2,c2)=mswap

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c      u1=1.0+q
c      u2=1.0+1.0/q
c      ddiv2=d/2.0

      subroutine 
     + invshoot(xst, yst, xit, yit, x1, x2, x3, y3, q1, q2, q3)
      implicit none
      real*8 xst, yst, xit, yit, x1, x2, x3, y3, q1, q2, q3
      real*8 r1, r2, r3 
      real*8 q1r1, q2r2, q3r3
       r1 = (xit-x1)**2+yit**2
       r2 = (xit-x2)**2+yit**2
       r3 = (xit-x3)**2+(yit-y3)**2
       q1r1 = q1/r1
       q2r2 = q2/r2
       q3r3 = q3/r3
       xst = xit - q1r1*(xit-x1) - 
     + q2r2*(xit-x2) -
     + q3r3*(xit-x3)
      yst = yit - q1r1*yit - 
     + q2r2*yit -
     + q3r3*(yit-y3)
       return
       end

      subroutine inverseshoot(xst,yst,xit,yit,ddiv2,u1,u2)

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c      u1=m1^(-1)
c      u2=m2^(-1)
       real*8 xst,yst,xit,yit,ddiv2,u1,u2
       real*8 xipd,ximd,xipd2,ximd2,yi2
       real*8 d1,d2,d3,d4
c      m1=1.0/(1.0+q)
c      m2=q/(1.0+q)
       xipd=xit+ddiv2
       ximd=xit-ddiv2
       xipd2=xipd*xipd
       ximd2=ximd*ximd
       yi2=yit*yit
       d1=xipd2+yi2
       d2=ximd2+yi2
       d3=d1*u1
       d4=d2*u2
       xst=xit-xipd/d3-ximd/d4
       yst=yit-yit/d3-yit/d4
c      xst=xit-m1*(xit+d/2.0)/((xit+d/2.)**2+yit**2)-
c     &m2*(xit-d/2.0)/((xit-d/2.)**2+yit**2)
c      yst=yit-m1*yit/((xit+d/2.)**2+yit**2)-
c     &m2*yit/((xit-d/2.)**2+yit**2)
      return
      end

      subroutine buildgrid(mini,maxi,origin,gridsize,
     +           gridmin,gridnum)
      implicit none
      real*8 mini,maxi,origin,gridsize,gridmin
      integer*4 gridnum,n1,n2

      n1=nint((mini-origin)/gridsize+0.5d0)
      n2=nint((maxi-origin)/gridsize-0.5d0)

      gridmin=n1*gridsize+origin
      gridnum=n2-n1+1

c      gridmin=origin+(int((mini-origin)/gridsize)-1)*gridsize
c      gridmax=origin+(int((maxi-origin)/gridsize)+1)*gridsize

c      gridnum=(gridmax-gridmin)/gridsize+1

      return
      end

      subroutine boundarytest(check,x1,y1,x2,y2,gstep,
     +                        rho2,scx,scy,
     + qx1, qx2, qx3, qy3, q1, q2, q3)

      implicit none
      integer*4 check,i
      real*8 x1,y1,x2,y2,x,y
      real*8 gstepsize
      integer*4 gstep
      real*8 rho2
      real*8 sxi,syi,scx,scy
      real*8 k1,k2
      real*8 qx1, qx2, qx3, qy3, q1, q2, q3

      if(gstep.LT.2) then
      check=1
      goto 40
      endif

      check=1
      k1=0d0
      k2=0d0

      if(x1.NE.x2) then
      k1=(y2-y1)/(x2-x1)
      endif

      if(y1.NE.y2) then
      k2=(x2-x1)/(y2-y1)
      endif

      if((x1.EQ.x2).AND.(y1.EQ.y2)) then
      check=1
      goto 40
      endif

      if (dabs(y1-y2).GT.dabs(x1-x2)) then

      gstepsize=(y2-y1)/gstep

      do 10 i=1,gstep-1

      y=y1+gstepsize*i
      x=x1+(y-y1)*k2

c      call inverseshoot(sxi,syi,x,y,bdiv2,u1,u2)
      call
     + invshoot(sxi, syi, x, y, qx1, qx2, qx3, qy3, q1, q2, q3)

      if (((sxi-scx)**2+(syi-scy)**2).LE.rho2) then
      check=0
      goto 40
      endif

10    continue

      else

      gstepsize=(x2-x1)/gstep

      do 20 i=1,gstep-1

      x=x1+gstepsize*i
      y=y1+(x-x1)*k1

c      call inverseshoot(sxi,syi,x,y,bdiv2,u1,u2)
      call
     + invshoot(sxi, syi, x, y, qx1, qx2, qx3, qy3, q1, q2, q3)

      if (((sxi-scx)**2+(syi-scy)**2).LT.rho2) then
      check=0
      goto 40
      endif

20    continue

      endif
40    return
      end


      subroutine crosstest(xa,ya,xb,yb,xc,yc,xd,yd,flag)
      implicit none
      real*8 xa,ya,xb,yb,xc,yc,xd,yd
      real*8 s1,s2,s3
      integer*4 flag
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     consider zero area                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      flag=1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     flag=1: not crossed                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      s1=xa*yb-ya*xb + xb*yc - yb*xc + xc*ya - yc*xa
      s2=xa*yb-ya*xb + xb*yd - yb*xd + xd*ya - yd*xa
      if(s1*s2.LT.0) then
      goto 10
      else
      s3=xa*yc-ya*xc + xc*yd - yc*xd + xd*ya - yd*xa
      if(s2*s3.LT.0) then
      goto 10
      else
      flag=0
      endif
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    check this to see the function of return            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
10    return
      end

      subroutine sort(array,dimen)
      implicit none
      integer*4 dimen
      real*8 array(dimen)
      integer*4 i,j
      real*8 swap

      do 10 i=2,dimen
         swap=array(i)
         do 20 j=i-1,1,-1
         if (array(j).LT.swap) then
         goto 30
         endif
      array(j+1)=array(j)
20    continue
      j=0
30    array(j+1)=swap
10    continue
      return
      end

      subroutine locate(xx,j1,n,x,j)
      INTEGER j,n
      real*8 x,xx(n)
      INTEGER j1,jm,ju
c      j1=0
      ju=n+1
10    if(ju-j1.gt.1) then
        jm=(ju+j1)/2
        if(x.ge.xx(jm)) then
           j1=jm
        else
           ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1)) then
        j=1
      else if(x.ge.xx(n)) then
        j=n
      else
        j=j1
      endif
      return
      END

c*********************************************************
c Numerical Recipe
c*********************************************************
      SUBROUTINE hpsorti(n,ra,rindex)
      integer*4 n
      real*8 ra(n)
      integer*4 rindex(n)
      real*8 i,ir,j,l
      real*8 rra
      integer*4 rrindex

      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
         if(l.gt.1)then
             l=l-1
             rra=ra(l)
             rrindex=rindex(l)
         else
             rra=ra(ir)
             ra(ir)=ra(1)

             rrindex=rindex(ir)
             rindex(ir)=rindex(1)

             ir=ir-1
             if(ir.eq.1)then
                 ra(1)=rra

                 rindex(1)=rrindex

                 return
             endif
          endif
             i=l
             j=l+l
20           if(j.le.ir)then
                if(j.lt.ir)then
                   if(ra(j).lt.ra(j+1))j=j+1
                endif
                if(rra.lt.ra(j))then
                   ra(i)=ra(j)

                    rindex(i)=rindex(j)

                    i=j
                    j=j+j
                else
                    j=ir+1
                endif
              goto 20
              endif
              ra(i)=rra

              rindex(i)=rrindex
      goto 10
      END



      subroutine reversearray (array1,column,n,dimen,max_sol)
      implicit none
      integer*4 max_sol
      integer*4 column,n,dimen,i
      real*8 array1(dimen,max_sol),temp1
 
      do 10 i=1,n/2
      
      temp1=array1(i,column)
      array1(i,column)=array1(n-i+1,column)
      array1(n-i+1,column)=temp1      

10    continue
      
      return
      end      
      
      subroutine stokes(rx,ry,lx,ly,n,col,area,dimen,max_sol)
      integer*4 max_sol
      integer*4 n,col,dimen,i
      real*8 rx(dimen*max_sol+max_sol,max_sol),
     + ry(dimen*max_sol+max_sol,max_sol)
      real*8 lx(dimen*max_sol+max_sol,max_sol),
     + ly(dimen*max_sol+max_sol,max_sol)
      real*8 area,p

      area=0
      do 10 i=1,n-1
      call cp(rx(i,col),lx(i,col),ry(i,col),ly(i,col),p)
      area=area+p/2
10    continue
      return
      end
     
      subroutine cp(x1,x2,y1,y2,p)
      real*8 x1,x2,y1,y2,p
      p=x1*y2-x2*y1
      return
      end
       
      subroutine dp(x1,x2,y1,y2,p)
      real*8 x1,x2,y1,y2,p
      p=x1*x2+y1*y2
      return
      end

      subroutine 
     + lenssolver(ixt,iyt,amt,msol,da,db,phi,ma,mb,sx,sy,magtott,
     + mparity,max_sol,errortype, dis_max,
     + x1, x2, x3, y3, q1, q2, q3,zr)
      implicit none
      real*8 da, db, phi, ma, mb, sx, sy, magtott
      integer max_sol, msol
      integer mparity(max_sol)
      real*8 ixt(max_sol), iyt(max_sol), amt(max_sol)
      integer errortype
      real*8 dis_max
      real*8 sx_temp, sy_temp, dis
      real*8 x1, x2, x3, y3, q1, q2, q3
      integer i
      integer nsolpos, nsolneg
      complex*16 zr(10)

      dis_max = -1.d0

      errortype = 0


      call getimt(ixt,iyt,amt,msol,da,db,phi,ma,mb,sx,sy,magtott,zr)

c      if (((msol.NE.10).AND.(msol.NE.8)
c     + .AND.(msol.NE.6).AND.(msol.NE.4))) then
c      errortype = 1
c      return
c      endif

      nsolpos=0
      nsolneg=0
      do i = 1,max_sol
      if(amt(i).gt.0) then
      mparity(i)=1
      nsolpos=nsolpos+1
      else if (amt(i).lt.0) then
      mparity(i)=-1
      nsolneg=nsolneg+1
      else
      mparity(i)=0
      endif
      enddo

c      if(((nsolneg-nsolpos).ne.2)) then
c      errortype = 1
c      return
c      endif

      do 30 i=1,max_sol
      if(mparity(i).ne.0) then
      call
     + invshoot
     + (sx_temp, sy_temp, ixt(i), iyt(i), x1, x2, x3, y3, q1, q2, q3)
      dis = sqrt((sx_temp - sx)**2 + (sy_temp - sy)**2)
      if(dis.gt.dis_max) then
      dis_max = dis
      endif
      endif
30    continue

      if(dis_max.lt.0) then
      errortype = 2
      return
      endif

      return
      end
