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
c     vx0=0.
c     call bdcnsx(0,margin,vx,vx0,ix)
      call bdfrex(0,margin,vx,ix)

      call bdfrex(1,margin,ro,ix)
c     vx0=0.
c     call bdcnsx(1,margin,vx,vx0,ix)
      call bdfrex(1,margin,vx,ix)


      return
      end
