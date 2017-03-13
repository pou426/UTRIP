c======================================================================|
      subroutine bnd(margin,ro,vx,vy,vz,bx,by,bz,az,ix,jx
     &           ,ipe,jpe,ipex,jpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx)
      dimension vx(ix,jx),vy(ix,jx),vz(ix,jx)
      dimension bx(ix,jx),by(ix,jx),bz(ix,jx),az(ix,jx)
c----------------------------------------------------------------------|      
      if (ipex.eq.1) then
      call bdperx(margin,margin,ro,ix,jx)
      call bdperx(margin,margin,vx,ix,jx)
      call bdperx(margin,margin,vy,ix,jx)
      call bdperx(margin,margin,vz,ix,jx)
      call bdperx(margin,margin,bx,ix,jx)
      call bdperx(margin,margin,by,ix,jx)
      call bdperx(margin,margin,bz,ix,jx)
      endif
      if (ipe.eq.0) then
      call bdfrex(0,margin,az,ix,jx)
      endif
      if (ipe.eq.ipex-1) then
      call bdfrex(1,margin,az,ix,jx)
      endif

      if (jpex.eq.1) then
      call bdpery(margin,margin,ro,ix,jx)
      call bdpery(margin,margin,vx,ix,jx)
      call bdpery(margin,margin,vy,ix,jx)
      call bdpery(margin,margin,vz,ix,jx)
      call bdpery(margin,margin,bx,ix,jx)
      call bdpery(margin,margin,by,ix,jx)
      call bdpery(margin,margin,bz,ix,jx)
      endif
      if (jpe.eq.0) then
      call bdfrey(0,margin,az,ix,jx)
      endif
      if (jpe.eq.jpex-1) then
      call bdfrey(1,margin,az,ix,jx)
      endif

      return
      end
