c======================================================================|
      subroutine bnd(margin,ro,vx,vy,by,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix),vx(ix),vy(ix),by(ix)
c----------------------------------------------------------------------|      
      call bdsppx(0,margin,ro,ix)
      call bdsmnx(0,margin,vx,ix)
      call bdsppx(0,margin,vy,ix)
      call bdsppx(0,margin,by,ix)

      call bdsppx(1,margin,ro,ix)
      call bdsmnx(1,margin,vx,ix)
      call bdsppx(1,margin,vy,ix)
      call bdsppx(1,margin,by,ix)


      return
      end
