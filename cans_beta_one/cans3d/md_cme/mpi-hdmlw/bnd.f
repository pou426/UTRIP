c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,vz,bx,by,bz,ix,jx,kx
     &           ,ipe,jpe,kpe,ipex,jpex,kpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
c----------------------------------------------------------------------|      

      if (ipe.eq.0) then
      call bdsppx(0,margin,ro,ix,jx,kx)
      call bdsppx(0,margin,pr,ix,jx,kx)
      call bdspnx(0,margin,vx,ix,jx,kx)
      call bdspnx(0,margin,vy,ix,jx,kx)
      call bdsppx(0,margin,vz,ix,jx,kx)
      call bdspnx(0,margin,bx,ix,jx,kx)
      call bdspnx(0,margin,by,ix,jx,kx)
      call bdsppx(0,margin,bz,ix,jx,kx)
      endif

      if (ipe.eq.ipex-1) then
      call bdfrex(1,margin,ro,ix,jx,kx)
      call bdfrex(1,margin,pr,ix,jx,kx)
      call bdfrex(1,margin,vx,ix,jx,kx)
      call bdfrex(1,margin,vy,ix,jx,kx)
      call bdfrex(1,margin,vz,ix,jx,kx)
      call bdfrex(1,margin,bx,ix,jx,kx)
      call bdfrex(1,margin,by,ix,jx,kx)
      call bdfrex(1,margin,bz,ix,jx,kx)
      endif

      if (jpe.eq.0) then
      call bdsppy(0,margin,ro,ix,jx,kx)
      call bdsppy(0,margin,pr,ix,jx,kx)
      call bdsppy(0,margin,vx,ix,jx,kx)
      call bdspny(0,margin,vy,ix,jx,kx)
      call bdspny(0,margin,vz,ix,jx,kx)
      call bdsppy(0,margin,bx,ix,jx,kx)
      call bdspny(0,margin,by,ix,jx,kx)
      call bdspny(0,margin,bz,ix,jx,kx)
      endif

      if (jpe.eq.jpex-1) then
      call bdsppy(1,margin,ro,ix,jx,kx)
      call bdsppy(1,margin,pr,ix,jx,kx)
      call bdsppy(1,margin,vx,ix,jx,kx)
      call bdspny(1,margin,vy,ix,jx,kx)
      call bdspny(1,margin,vz,ix,jx,kx)
      call bdsppy(1,margin,bx,ix,jx,kx)
      call bdspny(1,margin,by,ix,jx,kx)
      call bdspny(1,margin,bz,ix,jx,kx)
      endif

      if (kpex.eq.1) then
      call bdperz(margin,margin,ro,ix,jx,kx)
      call bdperz(margin,margin,pr,ix,jx,kx)
      call bdperz(margin,margin,vx,ix,jx,kx)
      call bdperz(margin,margin,vy,ix,jx,kx)
      call bdperz(margin,margin,vz,ix,jx,kx)
      call bdperz(margin,margin,bx,ix,jx,kx)
      call bdperz(margin,margin,by,ix,jx,kx)
      call bdperz(margin,margin,bz,ix,jx,kx)
      endif

      return
      end
