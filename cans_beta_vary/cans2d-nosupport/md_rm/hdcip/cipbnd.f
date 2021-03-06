c======================================================================|
      subroutine cipbnd(margin,ro,te,vxm,vym,rodx,tedx,vxdxm,vydxm
     &       ,rody,tedy,vxdym,vydym,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),te(ix,jx),vxm(ix,jx),vym(ix,jx)
      dimension rodx(ix,jx),tedx(ix,jx),vxdxm(ix,jx),vydxm(ix,jx)
      dimension rody(ix,jx),tedy(ix,jx),vxdym(ix,jx),vydym(ix,jx)
c----------------------------------------------------------------------|      
      call bdfrex(0,margin,ro,ix,jx)
      call bdfrex(0,margin,te,ix,jx)
      call bdfrex(0,margin-1,vxm,ix,jx)
      call bdfrex(0,margin-1,vym,ix,jx)
      call bdfrex(0,margin,rodx,ix,jx)
      call bdfrex(0,margin,tedx,ix,jx)
      call bdfrex(0,margin-1,vxdxm,ix,jx)
      call bdfrex(0,margin-1,vydxm,ix,jx)
      call bdfrex(0,margin,rody,ix,jx)
      call bdfrex(0,margin,tedy,ix,jx)
      call bdfrex(0,margin-1,vxdym,ix,jx)
      call bdfrex(0,margin-1,vydym,ix,jx)

      call bdfrex(1,margin,ro,ix,jx)
      call bdfrex(1,margin,te,ix,jx)
      call bdfrex(1,margin,vxm,ix,jx)
      call bdfrex(1,margin,vym,ix,jx)
      call bdfrex(1,margin,rodx,ix,jx)
      call bdfrex(1,margin,tedx,ix,jx)
      call bdfrex(1,margin,vxdxm,ix,jx)
      call bdfrex(1,margin,vydxm,ix,jx)
      call bdfrex(1,margin,rody,ix,jx)
      call bdfrex(1,margin,tedy,ix,jx)
      call bdfrex(1,margin,vxdym,ix,jx)
      call bdfrex(1,margin,vydym,ix,jx)

      call bdfrey(0,margin,ro,ix,jx)
      call bdfrey(0,margin,te,ix,jx)
      call bdfrey(0,margin-1,vxm,ix,jx)
      call bdfrey(0,margin-1,vym,ix,jx)
      call bdfrey(0,margin,rodx,ix,jx)
      call bdfrey(0,margin,tedx,ix,jx)
      call bdfrey(0,margin-1,vxdxm,ix,jx)
      call bdfrey(0,margin-1,vydxm,ix,jx)
      call bdfrey(0,margin,rody,ix,jx)
      call bdfrey(0,margin,tedy,ix,jx)
      call bdfrey(0,margin-1,vxdym,ix,jx)
      call bdfrey(0,margin-1,vydym,ix,jx)

      call bdfrey(1,margin,ro,ix,jx)
      call bdfrey(1,margin,te,ix,jx)
      call bdfrey(1,margin,vxm,ix,jx)
      call bdfrey(1,margin,vym,ix,jx)
      call bdfrey(1,margin,rodx,ix,jx)
      call bdfrey(1,margin,tedx,ix,jx)
      call bdfrey(1,margin,vxdxm,ix,jx)
      call bdfrey(1,margin,vydxm,ix,jx)
      call bdfrey(1,margin,rody,ix,jx)
      call bdfrey(1,margin,tedy,ix,jx)
      call bdfrey(1,margin,vxdym,ix,jx)
      call bdfrey(1,margin,vydym,ix,jx)

      return
      end
