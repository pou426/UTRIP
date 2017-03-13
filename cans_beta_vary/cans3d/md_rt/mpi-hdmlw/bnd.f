c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,vz,ix,jx,kx,dzm
     &       ,ipe,jpe,kpe,ipex,jpex,kpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx,kx)
      dimension pr(ix,jx,kx)
      dimension vx(ix,jx,kx)
      dimension vy(ix,jx,kx)
      dimension vz(ix,jx,kx)
      dimension dzm(kx)
c----------------------------------------------------------------------|      
      if (ipex .eq. 1) then
      call bdperx(1,margin,ro,ix,jx,kx)
      call bdperx(1,margin,pr,ix,jx,kx)
      call bdperx(1,margin,vx,ix,jx,kx)
      call bdperx(1,margin,vy,ix,jx,kx)
      call bdperx(1,margin,vz,ix,jx,kx)
      endif

      if (jpex .eq. 1) then
      call bdpery(1,margin,ro,ix,jx,kx)
      call bdpery(1,margin,pr,ix,jx,kx)
      call bdpery(1,margin,vx,ix,jx,kx)
      call bdpery(1,margin,vy,ix,jx,kx)
      call bdpery(1,margin,vz,ix,jx,kx)
      endif

      if (kpe.eq.0) then
      call bdfrez(0,margin,ro,ix,jx,kx)
      call bdfrdz(0,margin,pr,dzm,ix,jx,kx)
      call bdfrez(0,margin,vx,ix,jx,kx)
      call bdfrez(0,margin,vy,ix,jx,kx)
      call bdspnz(0,margin,vz,ix,jx,kx)
      endif

      if (kpe.eq.kpex-1) then
      call bdfrez(1,margin,ro,ix,jx,kx)
      call bdfrdz(1,margin,pr,dzm,ix,jx,kx)
      call bdfrez(1,margin,vx,ix,jx,kx)
      call bdfrez(1,margin,vy,ix,jx,kx)
      call bdspnz(1,margin,vz,ix,jx,kx)
      endif

      return
      end
