c======================================================================|
      subroutine engine(ro,pr,vr,vph,vz,br,bph,bz,r,rm
     &               ,gm,mu,dt,dr,drm,dph,dphm,dz,dzm,ix,jx,kx,mfdim)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension mfdim(3),mfdim3(3,3)
      dimension dr(ix),drm(ix)
      dimension dph(jx),dphm(jx)
      dimension dz(kx),dzm(kx)

      dimension ro(ix,jx,kx),pr(ix,jx,kx),eh(ix,jx,kx)
      dimension ee(ix,jx,kx),rr(ix,jx,kx),rph(ix,jx,kx),rz(ix,jx,kx)
      dimension vr(ix,jx,kx),vph(ix,jx,kx),vz(ix,jx,kx)
      dimension br(ix,jx,kx),bph(ix,jx,kx),bz(ix,jx,kx)
      dimension er(ix,jx,kx),eph(ix,jx,kx),ez(ix,jx,kx)

      dimension roh(ix,jx,kx),prh(ix,jx,kx),ehh(ix,jx,kx)
      dimension eeh(ix,jx,kx),rrh(ix,jx,kx),rphh(ix,jx,kx),rzh(ix,jx,kx)
      dimension vrh(ix,jx,kx),vphh(ix,jx,kx),vzh(ix,jx,kx)
      dimension brh(ix,jx,kx),bphh(ix,jx,kx),bzh(ix,jx,kx)
      dimension erh(ix,jx,kx),ephh(ix,jx,kx),ezh(ix,jx,kx)

      dimension ro2(ix,jx,kx),ee2(ix,jx,kx)
      dimension rr2(ix,jx,kx),rph2(ix,jx,kx),rz2(ix,jx,kx)
      dimension br2(ix,jx,kx),bph2(ix,jx,kx),bz2(ix,jx,kx)

      dimension r(ix),rm(ix)

      double precision mu
c----------------------------------------------------------------------|
c     mfdim=(/1,0,1/)
      call mfdimto3(mfdim,mfdim3)

      call vvtorv(rr,rph,rz,vr,vph,vz,ro,ix,jx,kx)
      call pr2eh(eh,pr,gm,ix,jx,kx)
      call eh2ee(ee,eh,ro,vr,vph,vz,br,bph,bz,mu,ix,jx,kx)

      call mlwreset(ro2,roh,ro,ix,jx,kx,mfdim)
      call mlwreset(ee2,eeh,ee,ix,jx,kx,mfdim)
      call mlwreset(rr2,rrh,rr,ix,jx,kx,mfdim)
      call mlwreset(rph2,rphh,rph,ix,jx,kx,mfdim)
      call mlwreset(rz2,rzh,rz,ix,jx,kx,mfdim)
      call mlwreset(br2,brh,br,ix,jx,kx,mfdim)
      call mlwreset(bph2,bphh,bph,ix,jx,kx,mfdim)
      call mlwreset(bz2,bzh,bz,ix,jx,kx,mfdim)

      call rvtovv(vr ,vph ,vz ,rr ,rph ,rz ,ro ,ix,jx,kx)
      call ee2eh(eh ,ee ,ro ,vr ,vph ,vz ,br ,bph ,bz ,mu,ix,jx,kx)
      call eh2pr(pr ,eh ,gm,ix,jx,kx)
      call ohmideal(er,eph,ez,br,bph,bz,vr,vph,vz,ix,jx,kx)
      call mlw_mh_c(ro ,pr ,eh ,vr ,vph ,vz ,br ,bph ,bz 
     &    ,er ,eph ,ez ,r
     &   ,ro2,ee2,rr2,rph2,rz2,br2,bph2,bz2
     &   ,roh,eeh,rrh,rphh,rzh,brh,bphh,bzh
     &   ,mu,dt,dr,drm,dph,dphm,dz,dzm,ix,jx,kx,mfdim,mfdim3,1)

      call rvtovv(vrh,vphh,vzh,rrh,rphh,rzh,roh,ix,jx,kx)
      call ee2eh(ehh,eeh,roh,vrh,vphh,vzh,brh,bphh,bzh,mu,ix,jx,kx)
      call eh2pr(prh,ehh,gm,ix,jx,kx)
      call ohmideal(erh,ephh,ezh,brh,bphh,bzh,vrh,vphh,vzh,ix,jx,kx)
      call mlw_mh_c(roh,prh,ehh,vrh,vphh,vzh,brh,bphh,bzh
     &    ,erh,ephh,ezh,rm
     &   ,ro2,ee2,rr2,rph2,rz2,br2,bph2,bz2
     &   ,roh,eeh,rrh,rphh,rzh,brh,bphh,bzh
     &   ,mu,dt,dr,drm,dph,dphm,dz,dzm,ix,jx,kx,mfdim,mfdim3,2)

      call rvtovv(vr,vph,vz,rr,rph,rz,ro,ix,jx,kx)
      qav=3.d0
      vvmin=1.d-4
      call avlap_8
     &    (ro2,ee2,rr2,rph2,rz2,br2,bph2,bz2
     &    ,ro ,ee ,rr ,rph ,rz ,br ,bph ,bz 
     &    ,vr,vph,vz,dt,qav,vvmin
     &    ,dr,drm,dph,dphm,dz,dzm,ix,jx,kx,mfdim)

      call daupdate(ro,ro2,ix,jx,kx)
      call daupdate(ee,ee2,ix,jx,kx)
      call daupdate(rr,rr2,ix,jx,kx)
      call daupdate(rph,rph2,ix,jx,kx)
      call daupdate(rz,rz2,ix,jx,kx)
      call daupdate(br,br2,ix,jx,kx)
      call daupdate(bph,bph2,ix,jx,kx)
      call daupdate(bz,bz2,ix,jx,kx)

      call rvtovv(vr,vph,vz,rr,rph,rz,ro,ix,jx,kx)
      call ee2eh(eh,ee,ro,vr,vph,vz,br,bph,bz,mu,ix,jx,kx)
      call eh2pr(pr,eh,gm,ix,jx,kx)

      return
      end
