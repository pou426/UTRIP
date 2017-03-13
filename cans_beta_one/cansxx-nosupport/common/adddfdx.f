      subroutine adddfdx(dan,f,dt,dx,dy,dz,ix,jx,kx,mdir)
c----------------------------------------------------------------------|
c     version 1.1 (2005/02/08 Yuji SATO)
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dan(ix,jx,kx)
      dimension f(ix,jx,kx)
      dimension dx(ix),dy(jx),dz(kx)

c----------------------------------------------------------------------|
c     define limiter functions
c     1. minmod limiter

c----------------------------------------------------------------------|

      if (mdir .eq. 1) then
        do k=1,kx
        do j=1,jx
        do i=2,ix-1
           dan(i,j,k)=dan(i,j,k)-dt*( (f(i,j,k)-f(i-1,j,k))/dx(i) )
        enddo
        enddo
        enddo
      endif

      if (mdir .eq. 2) then
        do k=1,kx
        do j=2,jx-1
        do i=1,ix
           dan(i,j,k)=dan(i,j,k)-dt*( (f(i,j,k)-f(i,j-1,k))/dy(j) )
        enddo
        enddo
        enddo
      endif

      if (mdir .eq. 3) then
        do k=2,kx-1
        do j=1,jx
        do i=1,ix
           dan(i,j,k)=dan(i,j,k)-dt*( (f(i,j,k)-f(i,j,k-1))/dz(k) )
        enddo
        enddo
        enddo
      endif

      return
      end
