c======================================================================|
      subroutine pertub(pr,x,ix)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix)
      dimension pr(ix)
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      wexp=0.3d0
      prexp=-0.1d0

      do i=1,ix
         pr(i) = pr(i)*(1+prexp*exp(-(x(i)/wexp)**2))
      enddo
      
      return
      end