c======================================================================|
      subroutine cipbnd(margin,te,vxm,vy,by,rodx,tedx,vxdxm,vydx
     &           ,vstar,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension te(ix),vxm(ix),vy(ix),by(ix)
      dimension rodx(ix),tedx(ix),vxdxm(ix),vydx(ix)
c----------------------------------------------------------------------|      
      zero=0.0d0
      te0=1.d0
      call bdcnsx(0,margin,rodx,zero,ix)
      call bdcnsx(0,margin,te,te0,ix)
      call bdcnsx(0,margin,tedx,zero,ix)
      call bdcnsx(0,margin,vy,vstar,ix)
      call bdcnsx(0,margin,vydx,zero,ix)
      call bdfrex(0,margin,by,ix)
c
c  Don't comment out !
c  If following lines are commented out, then,
c  artificial oscillation appears.
c
c     call bdfrex(0,margin-1,vxm,ix)
c     call bdcnsx(0,margin-1,vxdxm,zero,ix)

      call bdcnsx(1,margin-1,rodx,zero,ix)
      call bdfrex(1,margin,te,ix)
      call bdcnsx(1,margin-1,tedx,zero,ix)
      call bdfrex(1,margin,vxm,ix)
      call bdcnsx(1,margin-1,vxdxm,zero,ix)
      call bdfrex(1,margin,vy,ix)
      call bdcnsx(1,margin-1,vydx,zero,ix)
      call bdfrex(1,margin,by,ix)

      return
      end
