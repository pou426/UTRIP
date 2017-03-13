c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,ix,jx,ipe,jpe,ipex,jpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx)
c----------------------------------------------------------------------|

      if (ipe.eq.0) then
        call bdfrex(0,margin,ro,ix,jx)
        call bdfrex(0,margin,pr,ix,jx)
        call bdfrex(0,margin,vx,ix,jx)
        call bdfrex(0,margin,vy,ix,jx)
      endif

      if (ipe.eq.ipex-1) then
        call bdfrex(1,margin,ro,ix,jx)
        call bdfrex(1,margin,pr,ix,jx)
        call bdfrex(1,margin,vx,ix,jx)
        call bdfrex(1,margin,vy,ix,jx)
      endif

      if (jpe.eq.0) then
        call bdfrey(0,margin,ro,ix,jx)
        call bdfrey(0,margin,pr,ix,jx)
        call bdfrey(0,margin,vx,ix,jx)
        call bdfrey(0,margin,vy,ix,jx)
      endif

      if (jpe.eq.jpex-1) then
        call bdfrey(1,margin,ro,ix,jx)
        call bdfrey(1,margin,pr,ix,jx)
        call bdfrey(1,margin,vx,ix,jx)
        call bdfrey(1,margin,vy,ix,jx)
      endif

      return
      end

