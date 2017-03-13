c======================================================================|
      subroutine bnd(margin,ro,pr,vx,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix),pr(ix),vx(ix)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,ro,ix)
      call bdperx(margin,margin,pr,ix)
      call bdperx(margin,margin,vx,ix)

      return
      end
