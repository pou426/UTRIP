c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vz,bx,bz,ay,ix,jx
     &           ,ipe,jpe,ipex,jpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vz(ix,jx)
      dimension bx(ix,jx),bz(ix,jx)
      dimension ay(ix,jx)
c----------------------------------------------------------------------|      
      if (ipe.eq.0) then
      call bdsppx(0,margin,ro,ix,jx)
      call bdsppx(0,margin,pr,ix,jx)
      call bdspnx(0,margin,vx,ix,jx)
      call bdsppx(0,margin,vz,ix,jx)
      call bdspnx(0,margin,bx,ix,jx)
      call bdsppx(0,margin,bz,ix,jx)
      call bdsppx(0,margin,ay,ix,jx)
      endif

      if (ipe.eq.ipex-1) then
      call bdsppx(1,margin,ro,ix,jx)
      call bdsppx(1,margin,pr,ix,jx)
      call bdspnx(1,margin,vx,ix,jx)
      call bdsppx(1,margin,vz,ix,jx)
      call bdspnx(1,margin,bx,ix,jx)
      call bdsppx(1,margin,bz,ix,jx)
      call bdsppx(1,margin,ay,ix,jx)
      endif

      if (jpe.eq.0) then
      call bdsppy(0,margin,ro,ix,jx)
      call bdsppy(0,margin,pr,ix,jx)
      call bdsppy(0,margin,vx,ix,jx)
      call bdspny(0,margin,vz,ix,jx)
      call bdspny(0,margin,bx,ix,jx)
      call bdsppy(0,margin,bz,ix,jx)
      call bdsppy(0,margin,ay,ix,jx)
      endif

      if (jpe.eq.jpex-1) then
      call bdsppy(1,margin,ro,ix,jx)
      call bdsppy(1,margin,pr,ix,jx)
      call bdsppy(1,margin,vx,ix,jx)
      call bdspny(1,margin,vz,ix,jx)
      call bdspny(1,margin,bx,ix,jx)
      call bdsppy(1,margin,bz,ix,jx)
      call bdsppy(1,margin,ay,ix,jx)
      endif

      return
      end
