c======================================================================|
      subroutine roereset(da2,dah,da,ix,jx,kx)
c======================================================================|
c     numerical solver of mhd equations by roe method with muscl
c     for ideal 1d simulation (2nd order)
c     version 1.1 (2001/08/24 naoya fukuda)
c     version 1.2 (2005/02/08 Yuji SATO)
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)

      dimension da(ix,jx,kx),da2(ix,jx,kx),dah(ix,jx,kx)

c----------------------------------------------------------------------|
c     proceed half step
c     computation of 1st order flux f(i,l)
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
         dah(i,j,k)=da(i,j,k)
         da2(i,j,k)=da(i,j,k)
      enddo
      enddo
      enddo

      return
      end
