c======================================================================|
      subroutine cipadv(da,dadx,dady,dadz,u,v,w,isft,jsft,ksft
     &                    ,dt,dxm,dym,dzm,ix,jx,kx)
c======================================================================|
c
c NAME  cipadv
c
c PURPOSE
c    advance advective phase of CIP method
c
c INPUTS & OUTPUTS
c    da(ix,jx,kx): [double] physical variable
c    dadx(ix,jx,kx): [double] physical variable gradient
c    dady(ix,jx,kx): [double] physical variable gradient
c    dadz(ix,jx,kx): [double] physical variable gradient
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix,jx,kx) is the variable array defined at grid bounds
c
c    u(ix,jx,kx),v(ix,jx,kx),w(ix,jx,kx): [double] advection velocity
c    isft: [integer] 
c         0: if physical variable is defined at grid points
c         1: if physical variable is defined between grid points
c    dt: [double] delta time
c    dxm(ix),dym(jx),dzm(kx) : [double] grid spacing
c         NOTE: If physical variable is defined at grid points,
c               use dx(ix) instead of dxm(ix).
c    ix,jx,kx: [integer] dimension size
c
c HISTORY
c    written 2003-6-1 K. Takahashi based on T. Yokoyama's code
c
c----------------------------------------------------------------------|
      implicit real*8 (a-h,o-z)
      dimension dxm(ix),dym(jx),dzm(kx)
      dimension da(ix,jx,kx),dan(ix,jx,kx)
      dimension dadx(ix,jx,kx),dady(ix,jx,kx),dadz(ix,jx,kx)
      dimension u(ix,jx,kx),v(ix,jx,kx),w(ix,jx,kx)
      dimension dadxn(ix,jx,kx),dadyn(ix,jx,kx),dadzn(ix,jx,kx)
c----------------------------------------------------------------------|

      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
        if (u(i,j,k).le.0.) then
          iup=i+1
          dx1=+dxm(i+isft)
        else
          iup=i-1
          dx1=-dxm(i-1+isft)
        endif
        dx2=dx1**2
        dx3=dx1**3

        if (v(i,j,k).le.0.) then
          jup=j+1
          dy1=+dym(j+jsft)
        else
          jup=j-1
          dy1=-dym(j-1+jsft)
        endif
        dy2=dy1**2
        dy3=dy1**3

        if (w(i,j,k).le.0.) then
          kup=k+1
          dz1=+dzm(k+ksft)
        else
          kup=k-1
          dz1=-dzm(k-1+ksft)
        endif
        dz2=dz1**2
        dz3=dz1**3

        a01 =    (dadx(i,j,k)+dadx(iup,j,k))/dx2
     &      + 2.*(da  (i,j,k)-da  (iup,j,k))/dx3
        a11 = 3.*(da  (iup,j,k)-   da  (i,j,k))/dx2
     &          -(dadx(iup,j,k)+2.*dadx(i,j,k))/dx1
        a02 =    (dady(i,j,k)+dady(i,jup,k))/dy2
     &      + 2.*(da  (i,j,k)-da  (i,jup,k))/dy3
        a12 = 3.*(da  (i,jup,k)-   da  (i,j,k))/dy2
     &          -(dady(i,jup,k)+2.*dady(i,j,k))/dy1
        a03 =    (dadz(i,j,k)+dadz(i,j,kup))/dz2
     &      + 2.*(da  (i,j,k)-da  (i,j,kup))/dz3
        a13 = 3.*(da  (i,j,kup)-   da  (i,j,k))/dz2
     &          -(dadz(i,j,kup)+2.*dadz(i,j,k))/dz1

        b1 = da(i,j,k)-da(iup,j,k)-da(i,jup,k)+da(iup,jup,k)
        b2 = da(i,j,k)-da(i,jup,k)-da(i,j,kup)+da(i,jup,kup)
        b3 = da(i,j,k)-da(iup,j,k)-da(i,j,kup)+da(iup,j,kup)

        a05 = ( b1-(-dady(i,j,k)+dady(iup,j,k))*dy1)/(dx1*dy2)
        a04 = ( b1-(-dadx(i,j,k)+dadx(i,jup,k))*dx1)/(dx2*dy1)
        a14 = (-b1+(-dadx(i,j,k)+dadx(i,jup,k))*dx1
     &            +(-dady(i,j,k)+dady(iup,j,k))*dy1)/(dx1*dy1)

        a09 = ( b2-(-dadz(i,j,k)+dadz(i,jup,k))*dz1)/(dy1*dz2)
        a08 = ( b2-(-dady(i,j,k)+dady(i,j,kup))*dy1)/(dy2*dz1)
        a15 = (-b2+(-dady(i,j,k)+dady(i,j,kup))*dy1
     &            +(-dadz(i,j,k)+dadz(i,jup,k))*dz1)/(dy1*dz1)

        a06 = ( b3-(-dadz(i,j,k)+dadz(iup,j,k))*dz1)/(dx1*dz2)
        a07 = ( b3-(-dadx(i,j,k)+dadx(i,j,kup))*dx1)/(dx2*dz1)
        a16 = (-b3+(-dadz(i,j,k)+dadz(iup,j,k))*dz1
     &            +(-dadx(i,j,k)+dadx(i,j,kup))*dx1)/(dx1*dz1)

        a10 = (- da(i  ,j  ,k  )
     &         +(da(iup,j  ,k  )+da(i,jup,k  )+da(i  ,j,kup))
     &         -(da(iup,jup,k  )+da(i,jup,kup)+da(iup,j,kup))
     &         + da(iup,jup,kup))/(dx1*dy1*dz1)

        xx=-u(i,j,k)*dt
        yy=-v(i,j,k)*dt
        zz=-w(i,j,k)*dt

        dan(i,j,k)=((a01*xx+a04*yy+a07*zz+a11)*xx
     &            +a14*yy+dadx(i,j,k))*xx
     &            +((a02*yy+a05*xx+a08*zz+a12)*yy
     &            +a15*zz+dady(i,j,k))*yy
     &            +((a03*zz+a06*xx+a09*yy+a13)*zz
     &            +a16*xx+dadz(i,j,k))*zz
     &            +a10*xx*yy*zz+da(i,j,k)
        dadxn(i,j,k)=(3.*a01*xx+2.*(a04*yy+a07*zz+a11))*xx
     &              +(a05*yy+a10*zz+a14)*yy
     &              +(a06*zz+a16)*zz+dadx(i,j,k)
        dadyn(i,j,k)=(3.*a02*yy+2.*(a05*xx+a08*zz+a12))*yy
     &              +(a09*zz+a10*xx+a15)*zz
     &              +(a04*xx+a14)*xx+dady(i,j,k)
        dadzn(i,j,k)=(3.*a03*zz+2.*(a06*xx+a09*yy+a13))*zz
     &              +(a07*xx+a10*yy+a16)*xx
     &              +(a08*yy+a15)*yy+dadz(i,j,k)

        damin=min(da(i,j,k),da(iup,j,k),da(i,jup,k),da(i,j,kup))
        damax=max(da(i,j,k),da(iup,j,k),da(i,jup,k),da(i,j,kup))
        if (dan(i,j,k).lt.damin) dan(i,j,k)=damin
        if (dan(i,j,k).gt.damax) dan(i,j,k)=damax
      enddo
      enddo
      enddo

      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
        da(i,j,k)=dan(i,j,k)
        dadx(i,j,k)=dadxn(i,j,k)
        dady(i,j,k)=dadyn(i,j,k)
        dadz(i,j,k)=dadzn(i,j,k)
      enddo
      enddo
      enddo

      return
      end
