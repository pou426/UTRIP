c======================================================================|
      subroutine cipbnd(margin,vxm,vy,vz,by,bz
     &     ,rodx,vxdxm,vydx,vzdx,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension vxm(ix),vy(ix),vz(ix),by(ix),bz(ix)
      dimension rodx(ix),vxdxm(ix),vydx(ix),vzdx(ix)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,rodx,ix)
      call bdperx(margin-1,margin,vxm,ix)
      call bdperx(margin-1,margin,vxdxm,ix)
      call bdperx(margin,margin,vy,ix)
      call bdperx(margin,margin,vydx,ix)
      call bdperx(margin,margin,vz,ix)
      call bdperx(margin,margin,vzdx,ix)
      call bdperx(margin,margin,by,ix)
      call bdperx(margin,margin,bz,ix)


      return
      end
