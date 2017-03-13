c======================================================================|
      subroutine cipbnd(margin,ro,vxm,vym,rodx,vxdxm,vydxm
     &       ,rody,vxdym,vydym,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),vxm(ix,jx),vym(ix,jx)
      dimension rodx(ix,jx),vxdxm(ix,jx),vydxm(ix,jx)
      dimension rody(ix,jx),vxdym(ix,jx),vydym(ix,jx)
c----------------------------------------------------------------------|      
      call bdfrex(0,margin,ro,ix,jx)
      call bdcnsx(0,margin,rodx,0.d0,ix,jx)
      call bdfrex(0,margin,rody,ix,jx)
      call bdfrex(0,margin-1,vxm,ix,jx)
      call bdcnsx(0,margin-1,vxdxm,0.d0,ix,jx)
      call bdfrex(0,margin-1,vxdym,ix,jx)
      call bdfrex(0,margin,vym,ix,jx)
      call bdcnsx(0,margin,vydxm,0.d0,ix,jx)
      call bdfrex(0,margin,vydym,ix,jx)

      call bdfrex(1,margin,ro,ix,jx)
      call bdcnsx(1,margin,rodx,0.d0,ix,jx)
      call bdfrex(1,margin,rody,ix,jx)
      call bdfrex(1,margin,vxm,ix,jx)
      call bdcnsx(1,margin,vxdxm,0.d0,ix,jx)
      call bdfrex(1,margin,vxdym,ix,jx)
      call bdfrex(1,margin,vym,ix,jx)
      call bdcnsx(1,margin,vydxm,0.d0,ix,jx)
      call bdfrex(1,margin,vydym,ix,jx)

      call bdfrey(0,margin,ro,ix,jx)
      call bdfrey(0,margin,rodx,ix,jx)
      call bdcnsy(0,margin,rody,0.d0,ix,jx)
      call bdfrey(0,margin,vxm,ix,jx)
      call bdfrey(0,margin,vxdxm,ix,jx)
      call bdcnsy(0,margin,vxdym,0.d0,ix,jx)
      call bdfrey(0,margin-1,vym,ix,jx)
      call bdfrey(0,margin-1,vydxm,ix,jx)
      call bdcnsy(0,margin-1,vydym,0.d0,ix,jx)

      call bdfrey(1,margin,ro,ix,jx)
      call bdfrey(1,margin,rodx,ix,jx)
      call bdcnsy(1,margin,rody,0.d0,ix,jx)
      call bdfrey(1,margin,vxm,ix,jx)
      call bdfrey(1,margin,vxdxm,ix,jx)
      call bdcnsy(1,margin,vxdym,0.d0,ix,jx)
      call bdfrey(1,margin,vym,ix,jx)
      call bdfrey(1,margin,vydxm,ix,jx)
      call bdcnsy(1,margin,vydym,0.d0,ix,jx)

      return
      end
