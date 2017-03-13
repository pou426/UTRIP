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

      call bdsppx(1,margin,ro,ix,jx)
      call bdsppx(1,margin,pr,ix,jx)
      call bdspnx(1,margin,vx,ix,jx)
      call bdsppx(1,margin,vz,ix,jx)

      call bdsppy(0,margin,ro,ix,jx)
      call bdsppy(0,margin,pr,ix,jx)
      call bdsppy(0,margin,vx,ix,jx)
      call bdspny(0,margin,vz,ix,jx)

      call bdsppy(1,margin,ro,ix,jx)
      call bdsppy(1,margin,pr,ix,jx)
      call bdsppy(1,margin,vx,ix,jx)
      call bdspny(1,margin,vz,ix,jx)

      return
      end
