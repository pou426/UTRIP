c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,vz,bx,by,bz,droz,dprz,ix,jx,kx
     &           ,ipe,jpe,kpe,ipex,jpex,kpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension droz(ix,jx,margin),dprz(ix,jx,margin)
c----------------------------------------------------------------------|      
      if (ipex.eq.1) then
      call bdperx(margin,margin,pr,ix,jx,kx)
      call bdperx(margin,margin,ro,ix,jx,kx)
      call bdperx(margin,margin,vx,ix,jx,kx)
      call bdperx(margin,margin,vy,ix,jx,kx)
      call bdperx(margin,margin,vz,ix,jx,kx)
      call bdperx(margin,margin,bx,ix,jx,kx)
      call bdperx(margin,margin,by,ix,jx,kx)
      call bdperx(margin,margin,bz,ix,jx,kx)
      endif

      if (jpex.eq.1) then
      call bdpery(margin,margin,pr,ix,jx,kx)
      call bdpery(margin,margin,ro,ix,jx,kx)
      call bdpery(margin,margin,vx,ix,jx,kx)
      call bdpery(margin,margin,vy,ix,jx,kx)
      call bdpery(margin,margin,vz,ix,jx,kx)
      call bdpery(margin,margin,bx,ix,jx,kx)
      call bdpery(margin,margin,by,ix,jx,kx)
      call bdpery(margin,margin,bz,ix,jx,kx)
      endif

      if (kpe.eq.0) then
      call bdsppz(0,margin,pr,ix,jx,kx)
      call bdsppz(0,margin,ro,ix,jx,kx)
      call bdsppz(0,margin,vx,ix,jx,kx)
      call bdsppz(0,margin,vy,ix,jx,kx)
      call bdspnz(0,margin,vz,ix,jx,kx)
      call bdsppz(0,margin,bx,ix,jx,kx)
      call bdsppz(0,margin,by,ix,jx,kx)
      call bdspnz(0,margin,bz,ix,jx,kx)
      endif

      if (kpe.eq.kpex-1) then
      call bdcndz(1,margin,ro,droz,ix,jx,kx)
      call bdcndz(1,margin,pr,dprz,ix,jx,kx)
      call bdfrez(1,margin,vx,ix,jx,kx)
      call bdfrez(1,margin,vy,ix,jx,kx)
      call bdfrez(1,margin,vz,ix,jx,kx)
      call bdfrez(1,margin,bx,ix,jx,kx)
      call bdfrez(1,margin,by,ix,jx,kx)
      call bdfrez(1,margin,bz,ix,jx,kx)
      endif



      return
      end
