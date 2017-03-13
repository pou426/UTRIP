c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,vz,bx,by,bz,ix,jx,kx,mfdim)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension mfdim(3)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
c----------------------------------------------------------------------|      

      if (mfdim(1).eq.1) then
      call bdsppx(0,margin,ro,ix,jx,kx)
      call bdsppx(0,margin,pr,ix,jx,kx)
      call bdspnx(0,margin,vx,ix,jx,kx)
      call bdsppx(0,margin,vy,ix,jx,kx)
      call bdsppx(0,margin,vz,ix,jx,kx)
      call bdspnx(0,margin,bx,ix,jx,kx)
      call bdsppx(0,margin,by,ix,jx,kx)
      call bdsppx(0,margin,bz,ix,jx,kx)

      call bdsppx(1,margin,ro,ix,jx,kx)
      call bdsppx(1,margin,pr,ix,jx,kx)
      call bdspnx(1,margin,vx,ix,jx,kx)
      call bdsppx(1,margin,vy,ix,jx,kx)
      call bdsppx(1,margin,vz,ix,jx,kx)
      call bdspnx(1,margin,bx,ix,jx,kx)
      call bdsppx(1,margin,by,ix,jx,kx)
      call bdsppx(1,margin,bz,ix,jx,kx)
      endif

      if (mfdim(2).eq.1) then
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
      endif

      if (mfdim(3).eq.1) then
      call bdsppz(0,margin,ro,ix,jx,kx)
      call bdsppz(0,margin,pr,ix,jx,kx)
      call bdsppz(0,margin,vx,ix,jx,kx)
      call bdsppz(0,margin,vy,ix,jx,kx)
      call bdspnz(0,margin,vz,ix,jx,kx)
      call bdspnz(0,margin,bx,ix,jx,kx)
      call bdspnz(0,margin,by,ix,jx,kx)
      call bdsppz(0,margin,bz,ix,jx,kx)

      call bdsppz(1,margin,ro,ix,jx,kx)
      call bdsppz(1,margin,pr,ix,jx,kx)
      call bdsppz(1,margin,vx,ix,jx,kx)
      call bdsppz(1,margin,vy,ix,jx,kx)
      call bdspnz(1,margin,vz,ix,jx,kx)
      call bdspnz(1,margin,bx,ix,jx,kx)
      call bdspnz(1,margin,by,ix,jx,kx)
      call bdsppz(1,margin,bz,ix,jx,kx)
      endif

      return
      end
