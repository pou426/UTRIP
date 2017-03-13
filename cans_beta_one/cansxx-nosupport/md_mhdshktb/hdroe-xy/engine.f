c======================================================================|
      subroutine engine(ro,pr,vx,vy,vz,bx,by,bz
     &               ,gm,mu,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mfdim)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension mfdim(3),mfdim3(3,3)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension dz(kx),dzm(kx)

      dimension ro(ix,jx,kx),pr(ix,jx,kx),eh(ix,jx,kx)
      dimension ee(ix,jx,kx),rx(ix,jx,kx),ry(ix,jx,kx),rz(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension ex(ix,jx,kx),ey(ix,jx,kx),ez(ix,jx,kx)

      dimension roh(ix,jx,kx),prh(ix,jx,kx),ehh(ix,jx,kx)
      dimension eeh(ix,jx,kx),rxh(ix,jx,kx),ryh(ix,jx,kx),rzh(ix,jx,kx)
      dimension vxh(ix,jx,kx),vyh(ix,jx,kx),vzh(ix,jx,kx)
      dimension bxh(ix,jx,kx),byh(ix,jx,kx),bzh(ix,jx,kx)
      dimension exh(ix,jx,kx),eyh(ix,jx,kx),ezh(ix,jx,kx)

      dimension ro2(ix,jx,kx),ee2(ix,jx,kx)
      dimension rx2(ix,jx,kx),ry2(ix,jx,kx),rz2(ix,jx,kx)
      dimension bx2(ix,jx,kx),by2(ix,jx,kx),bz2(ix,jx,kx)

      double precision mu
c----------------------------------------------------------------------|

      pi = acos(-1.0d0)
      mu =4.d0*pi
c----------------------------------------------------------------------|
      call mfdimto3(mfdim,mfdim3)

      call vvtorv(rx,ry,rz,vx,vy,vz,ro,ix,jx,kx)
      call pr2eh(eh,pr,gm,ix,jx,kx)
      call eh2ee(ee,eh,ro,vx,vy,vz,bx,by,bz,mu,ix,jx,kx)

      call roereset(ro2,roh,ro,ix,jx,kx)
      call roereset(ee2,eeh,ee,ix,jx,kx)
      call roereset(rx2,rxh,rx,ix,jx,kx)
      call roereset(ry2,ryh,ry,ix,jx,kx)
      call roereset(rz2,rzh,rz,ix,jx,kx)
      call roereset(bx2,bxh,bx,ix,jx,kx)
      call roereset(by2,byh,by,ix,jx,kx)
      call roereset(bz2,bzh,bz,ix,jx,kx)

      dth=0.5d0*dt
      call roe_mh(ro ,pr ,vx ,vy ,vz ,bx ,by ,bz 
     &   ,roh,eeh,rxh,ryh,rzh,bxh,byh,bzh
     &   ,gm,dth,dx,dy,dz,ix,jx,kx,mfdim,1)

      call rvtovv(vxh,vyh,vzh,rxh,ryh,rzh,roh,ix,jx,kx)
      call ee2eh(ehh,eeh,roh,vxh,vyh,vzh,bxh,byh,bzh,mu,ix,jx,kx)
      call eh2pr(prh,ehh,gm,ix,jx,kx)

      call roe_mh(roh,prh,vxh,vyh,vzh,bxh,byh,bzh
     &   ,ro2,ee2,rx2,ry2,rz2,bx2,by2,bz2
     &   ,gm,dt,dx,dy,dz,ix,jx,kx,mfdim,2)

      call daupdate(ro,ro2,ix,jx,kx)
      call daupdate(ee,ee2,ix,jx,kx)
      call daupdate(rx,rx2,ix,jx,kx)
      call daupdate(ry,ry2,ix,jx,kx)
      call daupdate(rz,rz2,ix,jx,kx)
      call daupdate(bx,bx2,ix,jx,kx)
      call daupdate(by,by2,ix,jx,kx)
      call daupdate(bz,bz2,ix,jx,kx)

      call rvtovv(vx,vy,vz,rx,ry,rz,ro,ix,jx,kx)
      call ee2eh(eh,ee,ro,vx,vy,vz,bx,by,bz,mu,ix,jx,kx)
      call eh2pr(pr,eh,gm,ix,jx,kx)

      return
      end
