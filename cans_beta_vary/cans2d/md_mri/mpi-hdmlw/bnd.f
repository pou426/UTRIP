c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,vz,bx,by,bz,ay
     &   ,dvy,ix,jx
     &           ,ipe,jpe,ipex,jpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx),vz(ix,jx)
      dimension bx(ix,jx),by(ix,jx),bz(ix,jx),ay(ix,jx)
c----------------------------------------------------------------------|      
      if (ipex.eq.1) then
      call bdperx(margin,margin,pr,ix,jx)
      call bdperx(margin,margin,ro,ix,jx)
      call bdperx(margin,margin,vx,ix,jx)
      call bdperx(margin,margin,vy,ix,jx)
      call bdperx(margin,margin,vz,ix,jx)
      call bdperx(margin,margin,bx,ix,jx)
      call bdperx(margin,margin,by,ix,jx)
      call bdperx(margin,margin,bz,ix,jx)
      endif

      if (ipe.eq.0) then
      do i=1,margin
      do j=1,jx
        vy(i,j) =vy(i,j)-dvy
      enddo
      enddo
      call bdfrex(0,margin,ay,ix,jx)
      endif
      if (ipe.eq.ipex-1) then
      do i=1,margin
      do j=1,jx
        vy(ix-i+1,j) =vy(ix-i+1,j)+dvy
      enddo
      enddo
      call bdfrex(1,margin,ay,ix,jx)
      endif


      if (jpex.eq.1) then
      call bdpery(margin,margin,pr,ix,jx)
      call bdpery(margin,margin,ro,ix,jx)
      call bdpery(margin,margin,vx,ix,jx)
      call bdpery(margin,margin,vy,ix,jx)
      call bdpery(margin,margin,vz,ix,jx)
      call bdpery(margin,margin,bx,ix,jx)
      call bdpery(margin,margin,by,ix,jx)
      call bdpery(margin,margin,bz,ix,jx)
      endif
      if (jpe.eq.0) then
      call bdfrey(0,margin,ay,ix,jx)
      endif
      if (jpe.eq.jpex-1) then
      call bdfrey(1,margin,ay,ix,jx)
      endif

      return
      end
