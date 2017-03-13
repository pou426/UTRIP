c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vz,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx)
      dimension pr(ix,jx)
      dimension vx(ix,jx)
      dimension vz(ix,jx)
c----------------------------------------------------------------------|      
      call bdsppx(0,margin,ro,ix,jx)
      call bdsppx(0,margin,pr,ix,jx)
      call bdspnx(0,margin,vx,ix,jx)
      call bdsppx(0,margin,vz,ix,jx)

      call bdfrex(1,margin,ro,ix,jx)
      call bdfrex(1,margin,pr,ix,jx)
      call bdfrex(1,margin,vx,ix,jx)
      call bdfrex(1,margin,vz,ix,jx)

      call bdsmpy(0,margin,ro,ix,jx)
      call bdsmpy(0,margin,pr,ix,jx)
      call bdsmpy(0,margin,vx,ix,jx)
      call bdsmny(0,margin,vz,ix,jx)

      call bdfrey(1,margin,ro,ix,jx)
      call bdfrey(1,margin,pr,ix,jx)
      call bdfrey(1,margin,vx,ix,jx)
      call bdfrey(1,margin,vz,ix,jx)

      return
      end
