c======================================================================|
      subroutine bnd(margin,ro,vx,vy,vz,ix,jx,kx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx,kx)
      dimension vx(ix,jx,kx)
      dimension vy(ix,jx,kx)
      dimension vz(ix,jx,kx)
c----------------------------------------------------------------------|      
      call bdfrex(0,margin,ro,ix,jx,kx)
      call bdfrex(0,margin,vx,ix,jx,kx)
      call bdfrex(0,margin,vy,ix,jx,kx)
      call bdfrex(0,margin,vz,ix,jx,kx)

      call bdfrex(1,margin,ro,ix,jx,kx)
      call bdfrex(1,margin,vx,ix,jx,kx)
      call bdfrex(1,margin,vy,ix,jx,kx)
      call bdfrex(1,margin,vz,ix,jx,kx)

      call bdfrey(0,margin,ro,ix,jx,kx)
      call bdfrey(0,margin,vx,ix,jx,kx)
      call bdfrey(0,margin,vy,ix,jx,kx)
      call bdfrey(0,margin,vz,ix,jx,kx)

      call bdfrey(1,margin,ro,ix,jx,kx)
      call bdfrey(1,margin,vx,ix,jx,kx)
      call bdfrey(1,margin,vy,ix,jx,kx)
      call bdfrey(1,margin,vz,ix,jx,kx)

      call bdfrez(0,margin,ro,ix,jx,kx)
      call bdfrez(0,margin,vx,ix,jx,kx)
      call bdfrez(0,margin,vy,ix,jx,kx)
      call bdfrez(0,margin,vz,ix,jx,kx)

      call bdfrez(1,margin,ro,ix,jx,kx)
      call bdfrez(1,margin,vx,ix,jx,kx)
      call bdfrez(1,margin,vy,ix,jx,kx)
      call bdfrez(1,margin,vz,ix,jx,kx)

      return
      end

