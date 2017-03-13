c======================================================================|
       subroutine cipbnd(margin,ro,te,vxm,vym,bxm,bym
     &       ,rodx,tedx,vxdxm,vydxm
     &       ,rody,tedy,vxdym,vydym
     &       ,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),te(ix,jx),vxm(ix,jx),vym(ix,jx)
      dimension rodx(ix,jx),tedx(ix,jx),vxdxm(ix,jx),vydxm(ix,jx)
      dimension rody(ix,jx),tedy(ix,jx),vxdym(ix,jx),vydym(ix,jx)
      dimension bx(ix,jx),by(ix,jx),az(ix,jx)

      dimension tei(ix,jx),tedxi(ix,jx),tedyi(ix,jx)
      dimension rodxi(ix,jx),rodyi(ix,jx)

c----------------------------------------------------------------------|      
      call bdsppx(0,margin,ro,ix,jx)
      call bdsppx(0,margin,rodx,ix,jx)
      call bdspnx(0,margin,rody,ix,jx)
      call bdsppx(0,margin,te,ix,jx)
      call bdsppx(0,margin,tedx,ix,jx)
      call bdspnx(0,margin,tedy,ix,jx)
      call cipbdvsym(vxm,vxdxm,vxdym
     &              ,vym,vydxm,vydym
     &              ,margin,ix,jx,0,0,-1,+1)
      call cipbdbsym(bxm,bym
     &              ,margin,ix,jx,0,0,+1,-1)

      call bdsppx(1,margin,ro,ix,jx)
      call bdsppx(1,margin,rodx,ix,jx)
      call bdspnx(1,margin,rody,ix,jx)
      call bdsppx(1,margin,te,ix,jx)
      call bdsppx(1,margin,tedx,ix,jx)
      call bdspnx(1,margin,tedy,ix,jx)
      call cipbdvsym(vxm,vxdxm,vxdym
     &              ,vym,vydxm,vydym
     &              ,margin,ix,jx,0,1,-1,+1)
      call cipbdbsym(bxm,bym
     &              ,margin,ix,jx,0,1,-1,+1)

      call bdsppy(0,margin,ro,ix,jx)
      call bdsppy(0,margin,rodx,ix,jx)
      call bdspny(0,margin,rody,ix,jx)
      call bdsppy(0,margin,te,ix,jx)
      call bdsppy(0,margin,tedx,ix,jx)
      call bdspny(0,margin,tedy,ix,jx)
      call cipbdvsym(vxm,vxdxm,vxdym
     &              ,vym,vydxm,vydym
     &              ,margin,ix,jx,1,0,-1,+1)
      call cipbdbsym(bxm,bym
     &              ,margin,ix,jx,1,0,+1,-1)

      call bdsppy(1,margin,ro,ix,jx)
      call bdsppy(1,margin,rodx,ix,jx)
      call bdspny(1,margin,rody,ix,jx)
      call bdsppy(1,margin,te,ix,jx)
      call bdsppy(1,margin,tedx,ix,jx)
      call bdspny(1,margin,tedy,ix,jx)
      call cipbdvsym(vxm,vxdxm,vxdym
     &              ,vym,vydxm,vydym
     &              ,margin,ix,jx,1,1,-1,+1)
      call cipbdbsym(bxm,bym
     &              ,margin,ix,jx,1,1,+1,-1)

      return
      end
