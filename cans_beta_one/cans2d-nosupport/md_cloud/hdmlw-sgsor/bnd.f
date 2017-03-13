c======================================================================|
      subroutine bnd(margin,ro,vx,vy,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx)
      dimension vx(ix,jx)
      dimension vy(ix,jx)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,ro,ix,jx)
      call bdperx(margin,margin,vx,ix,jx)
      call bdperx(margin,margin,vy,ix,jx)

      call bdpery(margin,margin,ro,ix,jx)
      call bdpery(margin,margin,vx,ix,jx)
      call bdpery(margin,margin,vy,ix,jx)

      return
      end
