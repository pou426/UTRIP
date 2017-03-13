c======================================================================|
      subroutine glr_h(ro,pr,vx,x,dt,gm,dx,ix)
c======================================================================|
c
c NAME  glr_h
c
c PURPOSE
c    solve eqs. by Nonlinear Godunov type Lagrange Remap method with effects of
c        * hydrodynamics (van Leer 1979)
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    pr(ix): [double] pressure
c    vx(ix): [double] velocity 
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    dx(ix) : [double] grid spacing
c    gm: [double] polytropic index gamma
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-11-5 H. Koyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),x(ix)

      dimension ro(ix),pr(ix),vx(ix)
      dimension ee(ix),rx(ix),te(ix)
      dimension cs(ix),gradr(ix),gradp(ix),gradv(ix),pa(ix),vadt(ix)
      dimension xeuler(ix),dxinit(ix),dxinii(ix),dxlagc(ix)
      dimension advden(ix),advmox(ix),advene(ix)

      parameter( facdep = 0.5d0 ) 
      parameter( critif=1.0-1.0e-6 )
      gammp1=0.5*( gm + 1.0 )
      gammm1=0.5*( gm - 1.0 )
      rtgamm = sqrt( gm )
      gammx2=gammm1/rtgamm
      gammx1=gammm1/gm
c----------------------------------------------------------------------|
c     computation of conservative variables w(i,l)
      do i=1,ix
          cs(i) = sqrt(gm*pr(i)/ro(i))
          te(i) = pr(i)/(ro(i)*(gm-1.0d0))+0.5d0*vx(i)**2 
      enddo
      call monoto(dx,ro,gradr,ix)
      call monoto(dx,vx,gradv,ix)
      call monoto(dx,pr,gradp,ix)
c----------------------------------------------------------------------|
c     proceed Lagrange step
c----------------------------------------------------------------------|
      do i=2,ix
         depend=(dx(i-1)-cs(i-1)*dt)*facdep
         d1 = ro(i-1) + depend*gradr(i-1)
         p1 = pr(i-1) + depend*gradp(i-1)
         v1 = vx(i-1) + depend*gradv(i-1)
         depend=(dx(i)-cs(i)*dt)*facdep
         d2 = ro(i) - depend*gradr(i)
         p2 = pr(i) - depend*gradp(i)
         v2 = vx(i) - depend*gradv(i)

         rtpd1  = sqrt( p1*d1 )
         rtpd2  = sqrt( p2*d2 )
         paa    = (p1+p2)*0.5
         p1i    = 1.0/p1
         p2i    = 1.0/p2
         v1v2   = v1 - v2

c  iterate Rieman solver #1
         pasovp = paa*p2i
         phipop = sqrt( gammp1*pasovp + gammm1 )
         qmplus  = rtpd2*phipop

         pasovp = paa*p1i
         phipop = sqrt( gammp1*pasovp + gammm1 )
         qmminu = rtpd1*phipop

         paa     = ( qmplus*p1 + qmminu*p2 + v1v2*qmplus*qmminu )
     &           /( qmplus + qmminu )
         paa     = dmax1( paa, 0.d0 )

c  iterate Rieman solver #2
         pasovp = paa*p2i
         phipop = sqrt( gammp1*pasovp + gammm1 )
         qmplus  = rtpd2*phipop

         pasovp = paa*p1i
         phipop = sqrt( gammp1*pasovp + gammm1 )
         qmminu = rtpd1*phipop

         paa     = ( qmplus*p1 + qmminu*p2 + v1v2*qmplus*qmminu )
     &           /( qmplus + qmminu )
         paa     = dmax1( paa, 0.d0 )

c  iterate Rieman solver #3
         pasovp = paa*p2i
         phipop = sqrt( gammp1*pasovp + gammm1 )
         qmplus  = rtpd2*phipop

         pasovp = paa*p1i
         phipop = sqrt( gammp1*pasovp + gammm1 )
         qmminu = rtpd1*phipop

         paa     = ( qmplus*p1 + qmminu*p2 + v1v2*qmplus*qmminu )
     &           /( qmplus + qmminu )
         paa     = dmax1( paa, 0.d0 )

c  iterate Rieman solver #4
         pasovp = paa*p2i
         phipop = sqrt( gammp1*pasovp + gammm1 )
         qmplus  = rtpd2*phipop

         pasovp = paa*p1i
         phipop = sqrt( gammp1*pasovp + gammm1 )
         qmminu = rtpd1*phipop

         paa     = ( qmplus*p1 + qmminu*p2 + v1v2*qmplus*qmminu )
     &           /( qmplus + qmminu )
         paa     = dmax1( paa, 0.d0 )

         pa(i)   = paa
         vadt(i) = (qmminu*v1+qmplus*v2+p1-p2)/(qmplus+qmminu) * dt
      enddo
      
      do i=2,ix-1
         xeuler(i) = x(i)
         dxinit(i) = dx(i)
         dxinii(i) = 1.0d0/dx(i)
         dxlagc(i) = ro(i)*dx(i)
         dxlagi    = 1.d0/dxlagc(i)
         Vadt2 = vadt(i+1)
         Vadt1 = vadt(i)
         Pa2   = pa(i+1)
         Pa1   = pa(i)
c         write(*,*)i,pa(i),vadt(i),dt!pa1,pa2,vadt1,vadt2
         x(i)   = x(i) +0.5d0*(Vadt2 + Vadt1)
         vx(i)  = vx(i) -dt*dxlagi*(Pa2 - Pa1) 
         dx(i)  = dx(i) + (Vadt2 - Vadt1)
         volume = dx(i)*dxlagi
         ro(i)  = 1.0d0/volume
         te(i)  = te(i) - dxlagi*(Pa2*Vadt2 - Pa1*Vadt1)
c         speene = te(i)-0.5d0*vx(i)**2
c         pr(i)  = speene*(gm-1.0)*ro(i)
c         write(*,*)i,ro(i),te(i),vx(i)
      enddo
c      goto 999
c----------------------------------------------------------------------|
c     proceed Eulerian Remapping step
c----------------------------------------------------------------------|
      do i=2,ix-1
         rx(i) = ro(i)*vx(i)
         ee(i) = te(i)*ro(i)
      enddo

      call monoto(dx,ro,gradr,ix)
      call monoto(dx,rx,gradv,ix)
      call monoto(dx,ee,gradp,ix)

      do i=3,ix-1
         if(vadt(i).lt.0.0d0) then
            depend=(dx(i) + vadt(i))*facdep
            advden(i)=-vadt(i)*(ro(i)-depend*gradr(i))
            advmox(i)=-vadt(i)*(rx(i)-depend*gradv(i))
            advene(i)=-vadt(i)*(ee(i)-depend*gradp(i))
         else
            depend=(dx(i-1) - vadt(i))*facdep
            advden(i)=-vadt(i)*(ro(i-1)+depend*gradr(i-1))
            advmox(i)=-vadt(i)*(rx(i-1)+depend*gradv(i-1))
            advene(i)=-vadt(i)*(ee(i-1)+depend*gradp(i-1))
         endif
      enddo

      do i=3,ix-1
         ro(i) = (dxlagc(i)       + advden(i+1)-advden(i))*dxinii(i)
         rx(i) = (dxlagc(i)*vx(i) + advmox(i+1)-advmox(i))*dxinii(i) 
         ee(i) = (dxlagc(i)*te(i) + advene(i+1)-advene(i))*dxinii(i) 

         volume = 1.d0/ro(i)
         vx(i)  = rx(i)*volume
         x(i)   = xeuler(i)
         dx(i)  = dxinit(i)

         te(i)  = ee(i)*volume
         speene = te(i) - 0.5d0*vx(i)**2
         pr(i)  = ro(i)*(gm-1.0d0)*speene
      enddo
c----------------------------------------------------------------------|
 999  return
      end






