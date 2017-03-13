c======================================================================|
      subroutine bnd(margin,ro,pr,vx,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix)
      dimension pr(ix)
      dimension vx(ix)
c----------------------------------------------------------------------|      
c     call bdfrex(0,margin,ro,ix)
c     call bdfrex(0,margin,pr,ix)
c     call bdfrex(0,margin,vx,ix)

c     call bdfrex(1,margin,ro,ix)
c     call bdfrex(1,margin,pr,ix)
c     call bdfrex(1,margin,vx,ix)

      call bdsppx(0,margin,ro,ix)
      call bdsppx(0,margin,pr,ix)
      call bdspnx(0,margin,vx,ix)

      call bdsppx(1,margin,ro,ix)
      call bdsppx(1,margin,pr,ix)
      call bdspnx(1,margin,vx,ix)

      return
      end
