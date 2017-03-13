c======================================================================|
      subroutine bnd(margin,ro,pr,vx,ix,ipe,ipex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix)
      dimension pr(ix)
      dimension vx(ix)
c----------------------------------------------------------------------|      

      if (ipe.eq.0) then
      call bdsppx(0,margin,ro,ix)
      call bdsppx(0,margin,pr,ix)
      call bdspnx(0,margin,vx,ix)
      endif

      if (ipe.eq.ipex-1) then
      call bdsppx(1,margin,ro,ix)
      call bdsppx(1,margin,pr,ix)
      call bdspnx(1,margin,vx,ix)
      endif

      return
      end
