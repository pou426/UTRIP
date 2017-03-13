c======================================================================|
      subroutine bnd(margin,ey,ez,by,bz,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ey(ix),ez(ix),by(ix),bz(ix)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,ey,ix)
      call bdperx(margin,margin,ez,ix)
      call bdperx(margin,margin,by,ix)
      call bdperx(margin,margin,bz,ix)

      return
      end
