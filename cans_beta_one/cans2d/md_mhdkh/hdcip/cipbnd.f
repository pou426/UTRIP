c======================================================================|
      subroutine cipbnd(margin,ro,te,vxm,vym,vz,bxm,bym,bz
     &       ,rodx,tedx,vxdxm,vydxm,vzdx
     &       ,rody,tedy,vxdym,vydym,vzdy
     &       ,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),te(ix,jx),vxm(ix,jx),vym(ix,jx)
      dimension rodx(ix,jx),tedx(ix,jx),vxdxm(ix,jx),vydxm(ix,jx)
      dimension rody(ix,jx),tedy(ix,jx),vxdym(ix,jx),vydym(ix,jx)
      dimension bxm(ix,jx),bym(ix,jx),az(ix,jx)
      dimension vz(ix,jx),vzdx(ix,jx),vzdy(ix,jx)
      dimension bz(ix,jx)
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
      call bdperx(margin,margin,vz,ix,jx)
      call bdperx(margin,margin,vzdx,ix,jx)
      call bdperx(margin,margin,vzdy,ix,jx)
      call bdperx(margin,margin-1,bxm,ix,jx)
      call bdperx(margin,margin,bym,ix,jx)
      call bdperx(margin,margin,bz,ix,jx)

      call bdsppy(0,margin,ro,ix,jx)
      call bdsppy(0,margin,rodx,ix,jx)
      call bdspny(0,margin,rody,ix,jx)
      call bdsppy(0,margin,te,ix,jx)
      call bdsppy(0,margin,tedx,ix,jx)
      call bdspny(0,margin,tedy,ix,jx)
      call cipbdv3sym(vxm,vxdxm,vxdym
     &              ,vym,vydxm,vydym
     &              ,vz,vzdx,vzdy
     &              ,margin,ix,jx,1,0,-1,+1)
      call cipbdb3sym(bxm,bym,bz
     &              ,margin,ix,jx,1,0,+1,+1)

      call bdsppy(1,margin,ro,ix,jx)
      call bdsppy(1,margin,rodx,ix,jx)
      call bdspny(1,margin,rody,ix,jx)
      call bdsppy(1,margin,te,ix,jx)
      call bdsppy(1,margin,tedx,ix,jx)
      call bdspny(1,margin,tedy,ix,jx)
      call cipbdv3sym(vxm,vxdxm,vxdym
     &              ,vym,vydxm,vydym
     &              ,vz,vzdx,vzdy
     &              ,margin,ix,jx,1,1,-1,+1)
      call cipbdb3sym(bxm,bym,bz
     &              ,margin,ix,jx,1,1,+1,+1)

      return
      end
