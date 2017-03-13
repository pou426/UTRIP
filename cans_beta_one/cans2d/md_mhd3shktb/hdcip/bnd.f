c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,vz,bx,by,bz,az,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx),vz(ix,jx)
      dimension bx(ix,jx),by(ix,jx),bz(ix,jx),az(ix,jx)
c----------------------------------------------------------------------|      
      call bdfrex(0,margin,pr,ix,jx)
      call bdfrex(0,margin,ro,ix,jx)
      call bdfrex(0,margin,vx,ix,jx)
      call bdfrex(0,margin,vy,ix,jx)
      call bdfrex(0,margin,vz,ix,jx)
      call bdfrex(0,margin,bx,ix,jx)
      call bdfrex(0,margin,by,ix,jx)
      call bdfrex(0,margin,bz,ix,jx)
      call bdfrex(0,margin,az,ix,jx)

      call bdfrex(1,margin,pr,ix,jx)
      call bdfrex(1,margin,ro,ix,jx)
      call bdfrex(1,margin,vx,ix,jx)
      call bdfrex(1,margin,vy,ix,jx)
      call bdfrex(1,margin,vz,ix,jx)
      call bdfrex(1,margin,bx,ix,jx)
      call bdfrex(1,margin,by,ix,jx)
      call bdfrex(1,margin,bz,ix,jx)
      call bdfrex(1,margin,az,ix,jx)

      call bdfrey(0,margin,pr,ix,jx)
      call bdfrey(0,margin,ro,ix,jx)
      call bdfrey(0,margin,vx,ix,jx)
      call bdfrey(0,margin,vy,ix,jx)
      call bdfrey(0,margin,vz,ix,jx)
      call bdfrey(0,margin,bx,ix,jx)
      call bdfrey(0,margin,by,ix,jx)
      call bdfrey(0,margin,bz,ix,jx)
      call bdfrey(0,margin,az,ix,jx)

      call bdfrey(1,margin,pr,ix,jx)
      call bdfrey(1,margin,ro,ix,jx)
      call bdfrey(1,margin,vx,ix,jx)
      call bdfrey(1,margin,vy,ix,jx)
      call bdfrey(1,margin,vz,ix,jx)
      call bdfrey(1,margin,bx,ix,jx)
      call bdfrey(1,margin,by,ix,jx)
      call bdfrey(1,margin,bz,ix,jx)
      call bdfrey(1,margin,az,ix,jx)

      return
      end
