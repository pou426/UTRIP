c======================================================================|
      subroutine glr_h(ro,pr,vx,vy,x1,y1,dt,gm,dx1,ix,dy1,jx)
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
c    written 2002-11-6 H. Koyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx1(ix),x1(ix),dy1(jx),y1(jx)

      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension x(ix,jx),dx(ix,jx),v(ix,jx)
      
      dimension ee(ix,jx),rx(ix,jx),ry(ix,jx),te(ix,jx)
      dimension cs(ix,jx),gradr(ix,jx),gradp(ix,jx),gradv(ix,jx)
      dimension pa(ix,jx),vadt(ix,jx),gradu(ix,jx)
      dimension xeuler(ix,jx),dxinit(ix,jx),dxinii(ix,jx),dxlagc(ix,jx)
      dimension advden(ix,jx),advmox(ix,jx),advmoy(ix,jx),advene(ix,jx)

      parameter( facdep = 0.5d0 ) 
      parameter( critif=1.0-1.0e-6 )

c----------------------------------------------------------------------|
      gammp1=0.5*( gm + 1.0 )
      gammm1=0.5*( gm - 1.0 )
      rtgamm = sqrt( gm )
      gammx2=gammm1/rtgamm
      gammx1=gammm1/gm

      do m=1,2
        if (m.eq.1) then
          id=1
          jd=0
        else
          id=0
          jd=1
        endif

c----------------------------------------------------------------------|
c     computation of conservative variables w(i,l)

      if(id.eq.1)then
         do j=1,jx
         do i=1,ix
                x(i,j) =  x1(i)
               dx(i,j) = dx1(i)
                v(i,j) =  vx(i,j)
         enddo
         enddo
      endif

      if(jd.eq.1)then
         do j=1,jx
         do i=1,ix
            x(i,j) =  y1(j)
            dx(i,j) = dy1(j)
            v(i,j) =  vy(i,j)
         enddo
         enddo
      endif

      do j=1,jx
      do i=1,ix
          cs(i,j) = sqrt(gm*pr(i,j)/ro(i,j))
          te(i,j) = pr(i,j)/(ro(i,j)*(gm-1.0d0))
     &             +0.5d0*(vx(i,j)**2 + vy(i,j)**2) 
      enddo
      enddo

         call monoto(dx,ro,gradr,ix,jx,id,jd)
         call monoto(dx,v ,gradv,ix,jx,id,jd)
         call monoto(dx,pr,gradp,ix,jx,id,jd)
c----------------------------------------------------------------------|
c     proceed Lagrange step
c----------------------------------------------------------------------|
      do j=1+jd,jx
      do i=1+id,ix
         depend=(dx(i-id,j-jd)-cs(i-id,j-jd)*dt)*facdep
         d1 = ro(i-id,j-jd) + depend*gradr(i-id,j-jd)
         p1 = pr(i-id,j-jd) + depend*gradp(i-id,j-jd)
         v1 = v(i-id,j-jd)  + depend*gradv(i-id,j-jd)
         depend=(dx(i,j)-cs(i,j)*dt)*facdep
         d2 = ro(i,j) - depend*gradr(i,j)
         p2 = pr(i,j) - depend*gradp(i,j)
         v2 = v(i,j)  - depend*gradv(i,j)
c
         if(d1.le.0)write(*,*)'d1<=0'
         if(d2.le.0)write(*,*)'d2<=0'
         if(p1.le.0)write(*,*)'p1<=0',id,jd,i,j,p1,pr(i-id,j-jd)
         if(p2.le.0)write(*,*)'p2<=0',id,jd,i,j
c         
         rtpd1  = sqrt( p1*d1 )
         rtpd2  = sqrt( p2*d2 )
         paa    = (p1+p2)*0.5
         p1i    = 1.0/p1
         p2i    = 1.0/p2
         v1v2   = v1 - v2

c  iterate Rieman solver #1
         paold  = paa
         qmplus = rtpd2*sqrt( gammp1*paa*p2i + gammm1 )
         qmminu = rtpd1*sqrt( gammp1*paa*p1i + gammm1 )
         paa     = ( qmplus*p1 + qmminu*p2 + v1v2*qmplus*qmminu )
     &           /( qmplus + qmminu )
         paa     = dmax1( paa, paold*0.5d0 )

c  iterate Rieman solver #2
         paold  = paa
         qmplus = rtpd2*sqrt( gammp1*paa*p2i + gammm1 )
         qmminu = rtpd1*sqrt( gammp1*paa*p1i + gammm1 )
         paa     = ( qmplus*p1 + qmminu*p2 + v1v2*qmplus*qmminu )
     &           /( qmplus + qmminu )
         paa     = dmax1( paa, paold*0.5d0 )

c  iterate Rieman solver #3
         paold  = paa
         qmplus = rtpd2*sqrt( gammp1*paa*p2i + gammm1 )
         qmminu = rtpd1*sqrt( gammp1*paa*p1i + gammm1 )
         paa     = ( qmplus*p1 + qmminu*p2 + v1v2*qmplus*qmminu )
     &           /( qmplus + qmminu )
         paa     = dmax1( paa, paold*0.5d0 )

c  iterate Rieman solver #4
         paold  = paa
         qmplus = rtpd2*sqrt( gammp1*paa*p2i + gammm1 )
         qmminu = rtpd1*sqrt( gammp1*paa*p1i + gammm1 )
         paa     = ( qmplus*p1 + qmminu*p2 + v1v2*qmplus*qmminu )
     &           /( qmplus + qmminu )
         paa     = dmax1( paa, paold*0.5d0 )

         pa(i,j)   = paa
         vadt(i,j) = (qmminu*v1+qmplus*v2+p1-p2)/(qmplus+qmminu) * dt
         if(paa.le.0)write(*,*)'pa<=0',id,jd,i,j
      enddo
      enddo

      do j=1+jd,jx-jd
      do i=1+id,ix-id
         xeuler(i,j) = x(i,j)
         dxinit(i,j) = dx(i,j)
         dxinii(i,j) = 1.0d0/dxinit(i,j)
         dxlagc(i,j) = ro(i,j)*dxinit(i,j)
         dxlagi    = 1.d0/dxlagc(i,j)
         Vadt2 = vadt(i+id,j+jd)
         Vadt1 = vadt(i,j)
         Pa2   = pa(i+id,j+jd)
         Pa1   = pa(i,j)
c
         x(i,j)   = x(i,j) +0.5d0*(Vadt2 + Vadt1)
         v(i,j)   = v(i,j) -dt*dxlagi*(Pa2 - Pa1) 
         dx(i,j)  = dx(i,j) + (Vadt2 - Vadt1)
         volume   = dx(i,j)*dxlagi
         ro(i,j)  = 1.0d0/volume
         te(i,j)  = te(i,j) - dxlagi*(Pa2*Vadt2 - Pa1*Vadt1)
         if(te(i,j).le.0)write(*,*)'te<=0'
      enddo
      enddo

      if(id.eq.1)then
         do j=1+jd,jx-jd
         do i=1+id,ix-id
            vx(i,j) = v(i,j)
         enddo
         enddo
      endif
      if(jd.eq.1)then
         do j=1+jd,jx-jd
         do i=1+id,ix-id
            vy(i,j) = v(i,j)
         enddo
         enddo
      endif


         do j=1+jd,jx-jd
         do i=1+id,ix-id
            speene=te(i,j)-.5*(vx(i,j)**2+vy(i,j)**2)
            if(speene.le.0)write(*,*)'L:speene<=0'
         enddo
         enddo


c----------------------------------------------------------------------|
c     proceed Eulerian Remapping step
c----------------------------------------------------------------------|
      do j=1+jd,jx-jd
      do i=1+id,ix-id
         rx(i,j) = ro(i,j)*vx(i,j)
         ry(i,j) = ro(i,j)*vy(i,j)
         ee(i,j) = te(i,j)*ro(i,j)
      enddo
      enddo

         call monoto(dx,ro,gradr,ix,jx,id,jd)
         call monoto(dx,rx,gradv,ix,jx,id,jd)
         call monoto(dx,ry,gradu,ix,jx,id,jd)
         call monoto(dx,ee,gradp,ix,jx,id,jd)

      do j=1+2*jd,jx-jd
      do i=1+2*id,ix-id
         if(vadt(i,j).lt.0.0d0) then
            depend=(dx(i,j) + vadt(i,j))*facdep
            advden(i,j)=-vadt(i,j)*(ro(i,j)-depend*gradr(i,j))
            advmox(i,j)=-vadt(i,j)*(rx(i,j)-depend*gradv(i,j))
            advmoy(i,j)=-vadt(i,j)*(ry(i,j)-depend*gradu(i,j))
            advene(i,j)=-vadt(i,j)*(ee(i,j)-depend*gradp(i,j))
         else
         depend=(dx(i-id,j-jd) - vadt(i,j))*facdep
         advden(i,j)=-vadt(i,j)*(ro(i-id,j-jd)+depend*gradr(i-id,j-jd))
         advmox(i,j)=-vadt(i,j)*(rx(i-id,j-jd)+depend*gradv(i-id,j-jd))
         advmoy(i,j)=-vadt(i,j)*(ry(i-id,j-jd)+depend*gradu(i-id,j-jd))
         advene(i,j)=-vadt(i,j)*(ee(i-id,j-jd)+depend*gradp(i-id,j-jd))
         endif
      enddo
      enddo

      do j=1+2*jd,jx-2*jd
      do i=1+2*id,ix-2*id
         ro(i,j) = 
     & (dxlagc(i,j)        +advden(i+id,j+jd)-advden(i,j))*dxinii(i,j)
         rx(i,j) = 
     & (dxlagc(i,j)*vx(i,j)+advmox(i+id,j+jd)-advmox(i,j))*dxinii(i,j) 
         ry(i,j) = 
     & (dxlagc(i,j)*vy(i,j)+advmoy(i+id,j+jd)-advmoy(i,j))*dxinii(i,j) 
         ee(i,j) = 
     & (dxlagc(i,j)*te(i,j)+advene(i+id,j+jd)-advene(i,j))*dxinii(i,j) 
c
         if(ro(i,j).le.0)write(*,*)'ro<=0'
         if(ee(i,j).le.0)write(*,*)'ee<=0'
c
         volume = 1.d0/ro(i,j)
         vx(i,j)  = rx(i,j)*volume
         vy(i,j)  = ry(i,j)*volume
         x(i,j)   = xeuler(i,j)
         dx(i,j)  = dxinit(i,j)

         te(i,j)  = ee(i,j)*volume
         speene   = te(i,j) - 0.5d0*(vx(i,j)**2 + vy(i,j)**2)
         if(speene.le.0)write(*,*)'R:speene<=0',id,jd,i,j
         pr(i,j)  = ro(i,j)*(gm-1.0d0)*speene
      enddo
      enddo
c----------------------------------------------------------------------|
      enddo
      return
      end
