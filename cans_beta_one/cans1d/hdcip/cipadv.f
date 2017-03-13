c======================================================================|
      subroutine cipadv(da,dadx,u,isft,dt,dxm,ix)
c======================================================================|
c
c NAME  cipadv
c
c PURPOSE
c    advance advective phase of CIP method
c
c INPUTS & OUTPUTS
c    dadx(ix): [double] physical variable gradient
c    da(ix): [double] physical variable
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    u(ix): [double] advection velocity
c    isft: [integer] 
c         0: if physical variable is defined at grid points
c         1: if physical variable is defined between grid points
c    dt: [double] delta time
c    dxm(ix) : [double] grid spacing
c         NOTE: If physical variable is defined at grid points,
c               use dx(ix) instead of dxm(ix).
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension da(ix),dadx(ix),u(ix),dxm(ix)
      dimension dan(ix),dadxn(ix)
c----------------------------------------------------------------------|

      do i=2,ix-1

            xx=-u(i)*dt

          if (u(i).le.0.) then
            iup=i+1
            dx1=+dxm(i+isft)
          else
            iup=i-1
            dx1=-dxm(i-1+isft)
          endif

            dx2=dx1**2
            dx3=dx1**3

        a=(dadx(i)+dadx(iup))/dx2 + 2*(da(i)-da(iup))/dx3
        b=3.*(da(iup)-da(i))/dx2 - (2.*dadx(i)+dadx(iup))/dx1

        dan(i)=((a*xx + b)*xx + dadx(i))*xx + da(i)
        dadxn(i)=(3.*a*xx + 2.*b)*xx + dadx(i)

        damin=min(da(i),da(iup))
        damax=max(da(i),da(iup))
        if (dan(i).lt.damin) dan(i)=damin
        if (dan(i).gt.damax) dan(i)=damax
      enddo

      do i=1,ix
        da(i)=dan(i)
        dadx(i)=dadxn(i)
      enddo

      return
      end
