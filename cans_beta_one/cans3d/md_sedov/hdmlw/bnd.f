c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,vz,ix,jx,kx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
c----------------------------------------------------------------------|      

      call bdsppx(0,margin,ro,ix,jx,kx)
      call bdsppx(0,margin,pr,ix,jx,kx)
      call bdspnx(0,margin,vx,ix,jx,kx)
      call bdsppx(0,margin,vy,ix,jx,kx)
      call bdsppx(0,margin,vz,ix,jx,kx)

      call bdsppx(1,margin,ro,ix,jx,kx)
      call bdsppx(1,margin,pr,ix,jx,kx)
      call bdspnx(1,margin,vx,ix,jx,kx)
      call bdsppx(1,margin,vy,ix,jx,kx)
      call bdsppx(1,margin,vz,ix,jx,kx)

      call bdpery(margin,margin,ro,ix,jx,kx)
      call bdpery(margin,margin,pr,ix,jx,kx)
      call bdpery(margin,margin,vx,ix,jx,kx)
      call bdpery(margin,margin,vy,ix,jx,kx)
      call bdpery(margin,margin,vz,ix,jx,kx)

      call bdsmpz(0,margin,ro,ix,jx,kx)
      call bdsmpz(0,margin,pr,ix,jx,kx)
      call bdsmpz(0,margin,vx,ix,jx,kx)
      call bdsmpz(0,margin,vy,ix,jx,kx)
      call bdsmnz(0,margin,vz,ix,jx,kx)

      call bdsppz(1,margin,ro,ix,jx,kx)
      call bdsppz(1,margin,pr,ix,jx,kx)
      call bdsppz(1,margin,vx,ix,jx,kx)
      call bdsppz(1,margin,vy,ix,jx,kx)
      call bdspnz(1,margin,vz,ix,jx,kx)

      return
      end
