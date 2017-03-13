c======================================================================|
      subroutine bnd(margin,ro,vx,vy,bx,by,az,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx)
      dimension vx(ix,jx)
      dimension vy(ix,jx)
      dimension bx(ix,jx)
      dimension by(ix,jx)
      dimension az(ix,jx)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,ro,ix,jx)
      call bdperx(margin,margin,vx,ix,jx)
      call bdperx(margin,margin,vy,ix,jx)
      call bdperx(margin,margin,bx,ix,jx)
      call bdperx(margin,margin,by,ix,jx)
      call bdperx(margin,margin,az,ix,jx)

      call bdpery(margin,margin,ro,ix,jx)
      call bdpery(margin,margin,vx,ix,jx)
      call bdpery(margin,margin,vy,ix,jx)
      call bdpery(margin,margin,bx,ix,jx)
      call bdpery(margin,margin,by,ix,jx)
      call bdpery(margin,margin,az,ix,jx)

      return
      end
