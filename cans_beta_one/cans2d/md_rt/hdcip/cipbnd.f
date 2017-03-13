c======================================================================|
      subroutine cipbnd(margin,ro,te,vxm,vym,rodx,tedx,vxdxm,vydxm
     &       ,rody,tedy,vxdym,vydym,ix,jx,dym)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),te(ix,jx),vxm(ix,jx),vym(ix,jx)
      dimension rodx(ix,jx),tedx(ix,jx),vxdxm(ix,jx),vydxm(ix,jx)
      dimension rody(ix,jx),tedy(ix,jx),vxdym(ix,jx),vydym(ix,jx)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,ro,ix,jx)
      call bdperx(margin,margin,rodx,ix,jx)
      call bdperx(margin,margin,rody,ix,jx)
      call bdperx(margin,margin,te,ix,jx)
      call bdperx(margin,margin,tedx,ix,jx)
      call bdperx(margin,margin,tedy,ix,jx)
      call bdperx(margin,margin,vxm,ix,jx)
      call bdperx(margin,margin,vxdxm,ix,jx)
      call bdperx(margin,margin,vxdym,ix,jx)
      call bdperx(margin,margin,vym,ix,jx)
      call bdperx(margin,margin,vydxm,ix,jx)
      call bdperx(margin,margin,vydym,ix,jx)

      call bdfrey(0,margin,ro,ix,jx)
      call bdfrey(0,margin,rodx,ix,jx)
      call bdcnsy(0,margin,rody,0.d0,ix,jx)
      call bdfrdy(0,margin,te,dym,ix,jx)
      call bdfrey(0,margin,tedx,ix,jx)
      call bdfrey(0,margin,tedy,ix,jx)
      call bdfrey(0,margin,vxm,ix,jx)
      call bdfrey(0,margin,vxdxm,ix,jx)
      call bdcnsy(0,margin,vxdym,0.d0,ix,jx)
      call bdsmny(0,margin-1,vym,ix,jx)
      call bdsmpy(0,margin-1,vydxm,ix,jx)
      call bdcnsy(0,margin-1,vydym,0.d0,ix,jx)

      call bdfrey(1,margin,ro,ix,jx)
      call bdfrey(1,margin,rodx,ix,jx)
      call bdcnsy(1,margin,rody,0.d0,ix,jx)
      call bdfrdy(1,margin,te,dym,ix,jx)
      call bdfrey(1,margin,tedx,ix,jx)
      call bdfrey(1,margin,tedy,ix,jx)
      call bdfrey(1,margin,vxm,ix,jx)
      call bdfrey(1,margin,vxdxm,ix,jx)
      call bdcnsy(1,margin,vxdym,0.d0,ix,jx)
      call bdsmny(1,margin,vym,ix,jx)
      call bdsmpy(1,margin,vydxm,ix,jx)
      call bdcnsy(1,margin,vydym,0.d0,ix,jx)

      return
      end
