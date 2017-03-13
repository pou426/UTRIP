c======================================================================|
      subroutine bnd(margin,ro,vx,vy,vz,bx,by,bz,ix,jx,kx
     &       ,ipe,jpe,kpe,ipex,jpex,kpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
c----------------------------------------------------------------------|      
      if (ipex.eq.1) then
      call bdperx(margin,margin,ro,ix,jx,kx)
      call bdperx(margin,margin,vx,ix,jx,kx)
      call bdperx(margin,margin,vy,ix,jx,kx)
      call bdperx(margin,margin,vz,ix,jx,kx)
      call bdperx(margin,margin,bx,ix,jx,kx)
      call bdperx(margin,margin,by,ix,jx,kx)
      call bdperx(margin,margin,bz,ix,jx,kx)
      endif

      if (jpex.eq.1) then
      call bdpery(margin,margin,ro,ix,jx,kx)
      call bdpery(margin,margin,vx,ix,jx,kx)
      call bdpery(margin,margin,vy,ix,jx,kx)
      call bdpery(margin,margin,vz,ix,jx,kx)
      call bdpery(margin,margin,bx,ix,jx,kx)
      call bdpery(margin,margin,by,ix,jx,kx)
      call bdpery(margin,margin,bz,ix,jx,kx)
      endif

      if (kpex.eq.1) then
      call bdperz(margin,margin,ro,ix,jx,kx)
      call bdperz(margin,margin,vx,ix,jx,kx)
      call bdperz(margin,margin,vy,ix,jx,kx)
      call bdperz(margin,margin,vz,ix,jx,kx)
      call bdperz(margin,margin,bx,ix,jx,kx)
      call bdperz(margin,margin,by,ix,jx,kx)
      call bdperz(margin,margin,bz,ix,jx,kx)
      endif

      return
      end
