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
      call bdperx(margin,margin,ro,ix,jx)
      call bdperx(margin,margin,rodx,ix,jx)
      call bdperx(margin,margin,rody,ix,jx)
      call bdperx(margin,margin,te,ix,jx)
      call bdperx(margin,margin,tedx,ix,jx)
      call bdperx(margin,margin,tedy,ix,jx)
      call bdperx(margin,margin-1,vxm,ix,jx)
      call bdperx(margin,margin-1,vxdxm,ix,jx)
      call bdperx(margin,margin-1,vxdym,ix,jx)
      call bdperx(margin,margin,vym,ix,jx)
      call bdperx(margin,margin,vydxm,ix,jx)
      call bdperx(margin,margin,vydym,ix,jx)

      call bdsppy(0,margin,ro,ix,jx)
      call bdsppy(0,margin,rodx,ix,jx)
      call bdspny(0,margin,rody,ix,jx)
      call bdsppy(0,margin,te,ix,jx)
      call bdsppy(0,margin,tedx,ix,jx)
      call bdspny(0,margin,tedy,ix,jx)
      call cipbdvsym(vxm,vxdxm,vxdym
     &              ,vym,vydxm,vydym
     &              ,margin,ix,jx,1,0,-1,+1)

      call bdsppy(1,margin,ro,ix,jx)
      call bdsppy(1,margin,rodx,ix,jx)
      call bdspny(1,margin,rody,ix,jx)
      call bdsppy(1,margin,te,ix,jx)
      call bdsppy(1,margin,tedx,ix,jx)
      call bdspny(1,margin,tedy,ix,jx)
      call cipbdvsym(vxm,vxdxm,vxdym
     &              ,vym,vydxm,vydym
     &              ,margin,ix,jx,1,1,-1,+1)


      return
      end
