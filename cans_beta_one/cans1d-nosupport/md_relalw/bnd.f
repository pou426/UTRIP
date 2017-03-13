c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,by,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix),pr(ix),vx(ix),vy(ix),by(ix)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,pr,ix)
      call bdperx(margin,margin,ro,ix)
      call bdperx(margin,margin,vx,ix)
      call bdperx(margin,margin,vy,ix)
      call bdperx(margin,margin,by,ix)


      return
      end
