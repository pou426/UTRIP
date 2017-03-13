c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,vz,bx,by,bz,az,ix,jx
     &           ,ipe,jpe,ipex,jpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx),vz(ix,jx)
      dimension bx(ix,jx),by(ix,jx),bz(ix,jx)
      dimension az(ix,jx)
c----------------------------------------------------------------------|      
c     if (ipe.eq.0) then
c     call bdsppx(0,margin,ro,ix,jx)
c     call bdsppx(0,margin,pr,ix,jx)
c     call bdspnx(0,margin,vx,ix,jx)
c     call bdspnx(0,margin,vy,ix,jx)
c     call bdsppx(0,margin,vz,ix,jx)
c     call bdspnx(0,margin,bx,ix,jx)
c     call bdspnx(0,margin,by,ix,jx)
c     call bdsppx(0,margin,bz,ix,jx)
c     call bdsppx(0,margin,az,ix,jx)
c     endif

      if (ipe.eq.ipex-1) then
      call bdfrex(1,margin,ro,ix,jx)
      call bdfrex(1,margin,pr,ix,jx)
      call bdfrex(1,margin,vx,ix,jx)
      call bdfrex(1,margin,vy,ix,jx)
      call bdfrex(1,margin,vz,ix,jx)
      call bdfrex(1,margin,bx,ix,jx)
      call bdfrex(1,margin,by,ix,jx)
      call bdfrex(1,margin,bz,ix,jx)
      call bdfrex(1,margin,az,ix,jx)
      endif

      if (jpe.eq.0) then
      call bdsppy(0,margin,ro,ix,jx)
      call bdsppy(0,margin,pr,ix,jx)
      call bdsppy(0,margin,vx,ix,jx)
      call bdspny(0,margin,vy,ix,jx)
      call bdspny(0,margin,vz,ix,jx)
      call bdsppy(0,margin,bx,ix,jx)
      call bdspny(0,margin,by,ix,jx)
      call bdspny(0,margin,bz,ix,jx)
      call bdsppy(0,margin,az,ix,jx)
      endif

      if (jpe.eq.jpex-1) then
      call bdsppy(1,margin,ro,ix,jx)
      call bdsppy(1,margin,pr,ix,jx)
      call bdsppy(1,margin,vx,ix,jx)
      call bdspny(1,margin,vy,ix,jx)
      call bdspny(1,margin,vz,ix,jx)
      call bdsppy(1,margin,bx,ix,jx)
      call bdspny(1,margin,by,ix,jx)
      call bdspny(1,margin,bz,ix,jx)
      call bdsppy(1,margin,az,ix,jx)
      endif

      return
      end
