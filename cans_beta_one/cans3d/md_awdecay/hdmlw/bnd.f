c======================================================================|
      subroutine bnd(margin,ro,vx,vy,vz,bx,by,bz,ix,jx,kx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx,kx)
      dimension vx(ix,jx,kx)
      dimension vy(ix,jx,kx)
      dimension vz(ix,jx,kx)
      dimension bx(ix,jx,kx)
      dimension by(ix,jx,kx)
      dimension bz(ix,jx,kx)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,ro,ix,jx,kx)
      call bdperx(margin,margin,vx,ix,jx,kx)
      call bdperx(margin,margin,vy,ix,jx,kx)
      call bdperx(margin,margin,vz,ix,jx,kx)
      call bdperx(margin,margin,bx,ix,jx,kx)
      call bdperx(margin,margin,by,ix,jx,kx)
      call bdperx(margin,margin,bz,ix,jx,kx)

      call bdpery(margin,margin,ro,ix,jx,kx)
      call bdpery(margin,margin,vx,ix,jx,kx)
      call bdpery(margin,margin,vy,ix,jx,kx)
      call bdpery(margin,margin,vz,ix,jx,kx)
      call bdpery(margin,margin,bx,ix,jx,kx)
      call bdpery(margin,margin,by,ix,jx,kx)
      call bdpery(margin,margin,bz,ix,jx,kx)

      call bdperz(margin,margin,ro,ix,jx,kx)
      call bdperz(margin,margin,vx,ix,jx,kx)
      call bdperz(margin,margin,vy,ix,jx,kx)
      call bdperz(margin,margin,vz,ix,jx,kx)
      call bdperz(margin,margin,bx,ix,jx,kx)
      call bdperz(margin,margin,by,ix,jx,kx)
      call bdperz(margin,margin,bz,ix,jx,kx)

      return
      end
