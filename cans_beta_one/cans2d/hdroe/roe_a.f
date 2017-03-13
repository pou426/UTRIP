c======================================================================|
      subroutine roe_a(ro,dt,vx,vy,dx,ix,dy,jx)
c======================================================================|
c     numerical solver of mhd equations by roe method with muscl
c     for ideal 1d simulation (2nd order)
c     version 1.1 (2001/08/24 naoya fukuda)
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)

      dimension dx(ix),dy(jx)

      dimension ro(ix,jx),vx(ix,jx),vy(ix,jx)

      dimension fro(ix,jx)
      dimension roh(ix,jx)
      dimension row(ix,jx,2),vxw(ix,jx,2),vyw(ix,jx,2)

c----------------------------------------------------------------------|
c     proceed half step
c     computation of 1st order flux f(i,l)
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=ro(i,j)
      enddo
      enddo

c     x - direction
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix-1
         row(i,j,1)=ro(i,j)
         vxw(i,j,1)=vx(i,j)
         row(i,j,2)=ro(i+1,j)
         vxw(i,j,2)=vx(i+1,j)
      enddo
      enddo

      call roeflux_a(fro,row,vxw,ix,jx)

      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=roh(i,j)+0.5d0*dt*( (fro(i-1,j)-fro(i,j))/dx(i) )
      enddo
      enddo

c     y - direction
c----------------------------------------------------------------------|

      do j=1,jx-1
      do i=1,ix
         row(i,j,1)=ro(i,j)
         vyw(i,j,1)=vy(i,j)
         row(i,j,2)=ro(i,j+1)
         vyw(i,j,2)=vy(i,j+1)
      enddo
      enddo

      call roeflux_a(fro,row,vyw,ix,jx)

      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=roh(i,j)+0.5d0*dt*( (fro(i,j-1)-fro(i,j))/dy(j) )
      enddo
      enddo

c----------------------------------------------------------------------|
c     proceed full step
c     computation of 2nd order flux f(i,l)
c----------------------------------------------------------------------|
c     x - direction
c----------------------------------------------------------------------|
      call tvdminmod(1,roh,row,ix,jx)

      call roeflux_a(fro,row,vxw,ix,jx)

      do j=3,jx-2
      do i=3,ix-2
         ro(i,j)=ro(i,j)+dt*( (fro(i-1,j)-fro(i,j))/dx(i) )
      enddo
      enddo


c     y - direction
c----------------------------------------------------------------------|
      call tvdminmod(2,roh,row,ix,jx)

      call roeflux_a(fro,row,vyw,ix,jx)

      do j=3,jx-2
      do i=3,ix-2
         ro(i,j)=ro(i,j)+dt*( (fro(i,j-1)-fro(i,j))/dy(j) )
      enddo
      enddo
c----------------------------------------------------------------------|
      return
      end
