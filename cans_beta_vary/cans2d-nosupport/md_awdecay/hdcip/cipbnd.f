c======================================================================|
      subroutine cipbnd(margin,ro,vxm,vym,bxm,bym
     &       ,rodx,vxdxm,vydxm
     &       ,rody,vxdym,vydym
     &       ,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),vxm(ix,jx),vym(ix,jx)
      dimension rodx(ix,jx),vxdxm(ix,jx),vydxm(ix,jx)
      dimension rody(ix,jx),vxdym(ix,jx),vydym(ix,jx)
      dimension bx(ix,jx),by(ix,jx),az(ix,jx)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,ro,ix,jx)
      call bdperx(margin,margin,rodx,ix,jx)
      call bdperx(margin,margin,rody,ix,jx)
      call bdperx(margin,margin,vxm,ix,jx)
      call bdperx(margin,margin,vxdxm,ix,jx)
      call bdperx(margin,margin,vxdym,ix,jx)
      call bdperx(margin,margin,vym,ix,jx)
      call bdperx(margin,margin,vydxm,ix,jx)
      call bdperx(margin,margin,vydym,ix,jx)

      call bdpery(margin,margin,ro,ix,jx)
      call bdpery(margin,margin,rodx,ix,jx)
      call bdpery(margin,margin,rody,ix,jx)
      call bdpery(margin,margin,vxm,ix,jx)
      call bdpery(margin,margin,vxdxm,ix,jx)
      call bdpery(margin,margin,vxdym,ix,jx)
      call bdpery(margin,margin,vym,ix,jx)
      call bdpery(margin,margin,vydxm,ix,jx)
      call bdpery(margin,margin,vydym,ix,jx)


      return
      end
