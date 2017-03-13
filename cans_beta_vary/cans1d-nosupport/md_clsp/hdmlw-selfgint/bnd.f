c======================================================================|
      subroutine bnd(margin,ro,vx,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix)
      dimension vx(ix)
c----------------------------------------------------------------------|      
      call bdfrex(0,margin,ro,ix)
c     call bdcnsx(0,margin,vx,0.d0,ix)
      call bdspnx(0,margin,vx,ix)

      call bdfrex(1,margin,ro,ix)
      call bdfrex(1,margin,vx,ix)

c     call bdfrex(1,margin,ro,ix)
c     vx0=0.
c     call bdcnsx(1,margin,vx,vx0,ix)


      return
      end
