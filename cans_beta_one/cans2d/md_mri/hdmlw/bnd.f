c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,vz,bx,by,bz,ay
     &   ,dvy,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx),vz(ix,jx)
      dimension bx(ix,jx),by(ix,jx),bz(ix,jx),ay(ix,jx)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,pr,ix,jx)
      call bdperx(margin,margin,ro,ix,jx)
      call bdperx(margin,margin,vx,ix,jx)
      call bdperx(margin,margin,vy,ix,jx)
      call bdperx(margin,margin,vz,ix,jx)
      call bdperx(margin,margin,bx,ix,jx)
      call bdperx(margin,margin,by,ix,jx)
      call bdperx(margin,margin,bz,ix,jx)
      call bdfrex(0,margin,ay,ix,jx)
      call bdfrex(1,margin,ay,ix,jx)

      do i=1,margin
      do j=1,jx
        vy(i,j) =vy(i,j)-dvy
      enddo
      enddo
      do i=1,margin
      do j=1,jx
        vy(ix-i+1,j) =vy(ix-i+1,j)+dvy
      enddo
      enddo


      call bdpery(margin,margin,pr,ix,jx)
      call bdpery(margin,margin,ro,ix,jx)
      call bdpery(margin,margin,vx,ix,jx)
      call bdpery(margin,margin,vy,ix,jx)
      call bdpery(margin,margin,vz,ix,jx)
      call bdpery(margin,margin,bx,ix,jx)
      call bdpery(margin,margin,by,ix,jx)
      call bdpery(margin,margin,bz,ix,jx)
      call bdfrey(0,margin,ay,ix,jx)
      call bdfrey(1,margin,ay,ix,jx)

      return
      end
