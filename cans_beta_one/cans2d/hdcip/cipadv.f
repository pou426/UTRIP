c======================================================================|
      subroutine cipadv(da,dadx,dady,u,v,isft,jsft,dt,dxm,dym,ix,jx)
c======================================================================|
c
c NAME  cipadv
c
c PURPOSE
c    advance advective phase of CIP method
c
c INPUTS & OUTPUTS
c    da(ix,jx): [double] physical variable
c    dadx(ix,jx): [double] physical variable gradient
c    dady(ix,jx): [double] physical variable gradient
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix,jx) is the variable array defined at grid bounds
c
c    u(ix,jx),v(ix,jx): [double] advection velocity
c    isft: [integer] 
c         0: if physical variable is defined at grid points
c         1: if physical variable is defined between grid points
c    dt: [double] delta time
c    dxm(ix),dym(jx) : [double] grid spacing
c         NOTE: If physical variable is defined at grid points,
c               use dx(ix) instead of dxm(ix).
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dxm(ix),dym(jx)
      dimension da(ix,jx),dadx(ix,jx),dady(ix,jx),u(ix,jx),v(ix,jx)
      dimension dan(ix,jx),dadxn(ix,jx),dadyn(ix,jx)
c----------------------------------------------------------------------|

      do j=2,jx-1
      do i=2,ix-1

        if (u(i,j).le.0.) then
          iup=i+1
          dx1=+dxm(i+isft)
        else
          iup=i-1
          dx1=-dxm(i-1+isft)
        endif
        dx2=dx1**2
        dx3=dx1**3

        if (v(i,j).le.0.) then
          jup=j+1
          dy1=+dym(j+jsft)
        else
          jup=j-1
          dy1=-dym(j-1+jsft)
        endif
        dy2=dy1**2
        dy3=dy1**3

        a8=da(i,j)-da(iup,j)-da(i,jup)+da(iup,jup)
        g=(-a8+(dadx(i,jup)-dadx(i,j))*dx1+(dady(iup,j)-dady(i,j))*dy1)
     &         /(dx1*dy1)

        a=(dadx(i,j)+dadx(iup,j))/dx2 + 2*(da(i,j)-da(iup,j))/dx3
        c=(a8-(dadx(i,jup)-dadx(i,j))*dx1)/(dx2*dy1)
        e=3.*(da(iup,j)-da(i,j))/dx2 - (2.*dadx(i,j)+dadx(iup,j))/dx1

        b=(dady(i,j)+dady(i,jup))/dy2 + 2*(da(i,j)-da(i,jup))/dy3
        d=(a8-(dady(iup,j)-dady(i,j))*dy1)/(dx1*dy2)
        f=3.*(da(i,jup)-da(i,j))/dy2 - (2.*dady(i,j)+dady(i,jup))/dy1

        xx=-u(i,j)*dt
        yy=-v(i,j)*dt

        dan(i,j)= da(i,j) + g*xx*yy
     &          +((a*xx + c*yy + e)*xx + dadx(i,j))*xx
     &          +((b*yy + d*xx + f)*yy + dady(i,j))*yy
        dadxn(i,j)= dadx(i,j)
     &          +(3.*a*xx + 2.*c*yy + 2.*e)*xx
     &          +(g + d*yy)*yy
        dadyn(i,j)= dady(i,j)
     &          +(3.*b*yy + 2.*d*xx + 2.*f)*yy
     &          +(g + c*xx)*xx

        damin=min(da(i,j),da(iup,j),da(i,jup))
        damax=max(da(i,j),da(iup,j),da(i,jup))
        if (dan(i,j).lt.damin) dan(i,j)=damin
        if (dan(i,j).gt.damax) dan(i,j)=damax
      enddo
      enddo

      do j=2,jx-1
      do i=2,ix-1
        da(i,j)=dan(i,j)
        dadx(i,j)=dadxn(i,j)
        dady(i,j)=dadyn(i,j)
      enddo
      enddo

      return
      end
