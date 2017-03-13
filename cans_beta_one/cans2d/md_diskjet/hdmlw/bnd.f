c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,vz,bx,by,bz,ay,ix,jx,x,z)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx),vz(ix,jx)
      dimension bx(ix,jx),by(ix,jx),bz(ix,jx)
      dimension ay(ix,jx)
      dimension x(ix),z(jx)
c----------------------------------------------------------------------|      
      call bdsppx(0,margin,ro,ix,jx)
      call bdsppx(0,margin,pr,ix,jx)
      call bdspnx(0,margin,vx,ix,jx)
      call bdspnx(0,margin,vy,ix,jx)
      call bdsppx(0,margin,vz,ix,jx)
      call bdspnx(0,margin,bx,ix,jx)
      call bdspnx(0,margin,by,ix,jx)
      call bdsppx(0,margin,bz,ix,jx)
      call bdsppx(0,margin,ay,ix,jx)

      call bdfrex(1,margin,ro,ix,jx)
      call bdfrex(1,margin,pr,ix,jx)
      call bdfrex(1,margin,vx,ix,jx)
      call bdfrex(1,margin,vy,ix,jx)
      call bdfrex(1,margin,vz,ix,jx)
      call bdfrex(1,margin,bx,ix,jx)
      call bdfrex(1,margin,by,ix,jx)
      call bdfrex(1,margin,bz,ix,jx)
      call bdfrex(1,margin,ay,ix,jx)

      call bdsppy(0,margin,ro,ix,jx)
      call bdsppy(0,margin,pr,ix,jx)
      call bdsppy(0,margin,vx,ix,jx)
      call bdsppy(0,margin,vy,ix,jx)
      call bdspny(0,margin,vz,ix,jx)
      call bdspny(0,margin,bx,ix,jx)
      call bdspny(0,margin,by,ix,jx)
      call bdsppy(0,margin,bz,ix,jx)
      call bdsppy(0,margin,ay,ix,jx)

      call bdfrey(1,margin,ro,ix,jx)
      call bdfrey(1,margin,pr,ix,jx)
      call bdfrey(1,margin,vx,ix,jx)
      call bdfrey(1,margin,vy,ix,jx)
      call bdfrey(1,margin,vz,ix,jx)
      call bdfrey(1,margin,bx,ix,jx)
      call bdfrey(1,margin,by,ix,jx)
      call bdfrey(1,margin,bz,ix,jx)
      call bdfrey(1,margin,ay,ix,jx)


      return
      end
