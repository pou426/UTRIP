c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,vz,bx,by,bz,ix,jx,kx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
c----------------------------------------------------------------------|      
      call bdfrex(0,margin,ro,ix,jx,kx)
      call bdfrex(0,margin,pr,ix,jx,kx)
      call bdfrex(0,margin,vx,ix,jx,kx)
      call bdfrex(0,margin,vy,ix,jx,kx)
      call bdfrex(0,margin,vz,ix,jx,kx)
      call bdfrex(0,margin,bx,ix,jx,kx)
      call bdfrex(0,margin,by,ix,jx,kx)
      call bdfrex(0,margin,bz,ix,jx,kx)

      call bdfrex(1,margin,ro,ix,jx,kx)
      call bdfrex(1,margin,pr,ix,jx,kx)
      call bdfrex(1,margin,vx,ix,jx,kx)
      call bdfrex(1,margin,vy,ix,jx,kx)
      call bdfrex(1,margin,vz,ix,jx,kx)
      call bdfrex(1,margin,bx,ix,jx,kx)
      call bdfrex(1,margin,by,ix,jx,kx)
      call bdfrex(1,margin,bz,ix,jx,kx)

      call bdfrey(0,margin,ro,ix,jx,kx)
      call bdfrey(0,margin,pr,ix,jx,kx)
      call bdfrey(0,margin,vx,ix,jx,kx)
      call bdfrey(0,margin,vy,ix,jx,kx)
      call bdfrey(0,margin,vz,ix,jx,kx)
      call bdfrey(0,margin,bx,ix,jx,kx)
      call bdfrey(0,margin,by,ix,jx,kx)
      call bdfrey(0,margin,bz,ix,jx,kx)

      call bdfrey(1,margin,ro,ix,jx,kx)
      call bdfrey(1,margin,pr,ix,jx,kx)
      call bdfrey(1,margin,vx,ix,jx,kx)
      call bdfrey(1,margin,vy,ix,jx,kx)
      call bdfrey(1,margin,vz,ix,jx,kx)
      call bdfrey(1,margin,bx,ix,jx,kx)
      call bdfrey(1,margin,by,ix,jx,kx)
      call bdfrey(1,margin,bz,ix,jx,kx)

      call bdfrez(0,margin,ro,ix,jx,kx)
      call bdfrez(0,margin,pr,ix,jx,kx)
      call bdfrez(0,margin,vx,ix,jx,kx)
      call bdfrez(0,margin,vy,ix,jx,kx)
      call bdfrez(0,margin,vz,ix,jx,kx)
      call bdfrez(0,margin,bx,ix,jx,kx)
      call bdfrez(0,margin,by,ix,jx,kx)
      call bdfrez(0,margin,bz,ix,jx,kx)

      call bdfrez(1,margin,ro,ix,jx,kx)
      call bdfrez(1,margin,pr,ix,jx,kx)
      call bdfrez(1,margin,vx,ix,jx,kx)
      call bdfrez(1,margin,vy,ix,jx,kx)
      call bdfrez(1,margin,vz,ix,jx,kx)
      call bdfrez(1,margin,bx,ix,jx,kx)
      call bdfrez(1,margin,by,ix,jx,kx)
      call bdfrez(1,margin,bz,ix,jx,kx)

      return
      end
