c======================================================================|
      subroutine bnd(margin,ro,vx,vy,by,vz,bz,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix), vx(ix), vy(ix),by(ix), vz(ix),bz(ix)
c----------------------------------------------------------------------|      

      call bdperx(margin,margin,ro,ix)
      call bdperx(margin,margin,vx,ix)
      call bdperx(margin,margin,vy,ix)
      call bdperx(margin,margin,by,ix)
      call bdperx(margin,margin,vz,ix)
      call bdperx(margin,margin,bz,ix)

      return
      end
