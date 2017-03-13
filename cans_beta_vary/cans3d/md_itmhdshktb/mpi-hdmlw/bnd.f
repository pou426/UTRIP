c======================================================================|
      subroutine bnd(margin,ro,vx,vy,vz,bx,by,bz,ix,jx,kx
     &           ,ipe,jpe,kpe,ipex,jpex,kpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
c----------------------------------------------------------------------|      
      if (ipe.eq.0) then
      call bdfrex(0,margin,ro,ix,jx,kx)
      call bdfrex(0,margin,vx,ix,jx,kx)
      call bdfrex(0,margin,vy,ix,jx,kx)
      call bdfrex(0,margin,vz,ix,jx,kx)
      call bdfrex(0,margin,bx,ix,jx,kx)
      call bdfrex(0,margin,by,ix,jx,kx)
      call bdfrex(0,margin,bz,ix,jx,kx)
      endif

      if (ipe.eq.ipex-1) then
      call bdfrex(1,margin,ro,ix,jx,kx)
      call bdfrex(1,margin,vx,ix,jx,kx)
      call bdfrex(1,margin,vy,ix,jx,kx)
      call bdfrex(1,margin,vz,ix,jx,kx)
      call bdfrex(1,margin,bx,ix,jx,kx)
      call bdfrex(1,margin,by,ix,jx,kx)
      call bdfrex(1,margin,bz,ix,jx,kx)
      endif

      if (jpe.eq.0) then
      call bdfrey(0,margin,ro,ix,jx,kx)
      call bdfrey(0,margin,vx,ix,jx,kx)
      call bdfrey(0,margin,vy,ix,jx,kx)
      call bdfrey(0,margin,vz,ix,jx,kx)
      call bdfrey(0,margin,bx,ix,jx,kx)
      call bdfrey(0,margin,by,ix,jx,kx)
      call bdfrey(0,margin,bz,ix,jx,kx)
      endif

      if (jpe.eq.jpex-1) then
      call bdfrey(1,margin,ro,ix,jx,kx)
      call bdfrey(1,margin,vx,ix,jx,kx)
      call bdfrey(1,margin,vy,ix,jx,kx)
      call bdfrey(1,margin,vz,ix,jx,kx)
      call bdfrey(1,margin,bx,ix,jx,kx)
      call bdfrey(1,margin,by,ix,jx,kx)
      call bdfrey(1,margin,bz,ix,jx,kx)
      endif

      if (kpe.eq.0) then
      call bdfrez(0,margin,ro,ix,jx,kx)
      call bdfrez(0,margin,vx,ix,jx,kx)
      call bdfrez(0,margin,vy,ix,jx,kx)
      call bdfrez(0,margin,vz,ix,jx,kx)
      call bdfrez(0,margin,bx,ix,jx,kx)
      call bdfrez(0,margin,by,ix,jx,kx)
      call bdfrez(0,margin,bz,ix,jx,kx)
      endif

      if (kpe.eq.kpex-1) then
      call bdfrez(1,margin,ro,ix,jx,kx)
      call bdfrez(1,margin,vx,ix,jx,kx)
      call bdfrez(1,margin,vy,ix,jx,kx)
      call bdfrez(1,margin,vz,ix,jx,kx)
      call bdfrez(1,margin,bx,ix,jx,kx)
      call bdfrez(1,margin,by,ix,jx,kx)
      call bdfrez(1,margin,bz,ix,jx,kx)
      endif

      return
      end