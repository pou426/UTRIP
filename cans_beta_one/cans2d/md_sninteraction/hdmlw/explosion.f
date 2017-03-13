c======================================================================|
      subroutine explosion(pr,x,z,zexp,ix,jx)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),z(jx)
      dimension pr(ix,jx)

c----------------------------------------------------------------------|
      prism=1.e-8
      wexp=0.1d0

      do j=1,jx
      do i=1,ix
         ss=sqrt(x(i)**2+(z(j)-zexp)**2)
         pr(i,j)  = pr(i,j)+exp(-(ss/wexp)**2)
      enddo
      enddo

      
      return
      end
