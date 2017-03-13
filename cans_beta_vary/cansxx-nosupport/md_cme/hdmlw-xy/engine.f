c======================================================================|
      subroutine engine(ro,pr,vx,vy,vz,bx,by,bz,x,xm,y,ym,gx,gxm
     &               ,gm,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mfdim)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension mfdim(3),mfdim3(3,3)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension dz(kx),dzm(kx)
      dimension x(ix),xm(ix),y(jx),ym(jx)

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

      dimension gx(ix,jx,kx),gxm(ix,jx,kx)

      double precision mu
c----------------------------------------------------------------------|

      pi = acos(-1.0d0)
      mu =4.d0*pi
c----------------------------------------------------------------------|
      call mfdimto3(mfdim,mfdim3)

      call vvtorv(rx,ry,rz,vx,vy,vz,ro,ix,jx,kx)
      call pr2eh(eh,pr,gm,ix,jx,kx)
      call eh2ee(ee,eh,ro,vx,vy,vz,bx,by,bz,mu,ix,jx,kx)

      call mlwreset(ro2,roh,ro,ix,jx,kx,mfdim)
      call mlwreset(ee2,eeh,ee,ix,jx,kx,mfdim)
      call mlwreset(rx2,rxh,rx,ix,jx,kx,mfdim)
      call mlwreset(ry2,ryh,ry,ix,jx,kx,mfdim)
      call mlwreset(rz2,rzh,rz,ix,jx,kx,mfdim)
      call mlwreset(bx2,bxh,bx,ix,jx,kx,mfdim)
      call mlwreset(by2,byh,by,ix,jx,kx,mfdim)
      call mlwreset(bz2,bzh,bz,ix,jx,kx,mfdim)

      call rvtovv(vx ,vy ,vz ,rx ,ry ,rz ,ro ,ix,jx,kx)
      call ee2eh(eh ,ee ,ro ,vx ,vy ,vz ,bx ,by ,bz ,mu,ix,jx,kx)
      call eh2pr(pr ,eh ,gm,ix,jx,kx)
      call ohmideal(ex,ey,ez,bx,by,bz,vx,vy,vz,ix,jx,kx)
      call mlw_mh_s(ro ,pr ,eh ,vx ,vy ,vz ,bx ,by ,bz ,ex ,ey ,ez 
     &   ,x,y
     &   ,ro2,ee2,rx2,ry2,rz2,bx2,by2,bz2
     &   ,roh,eeh,rxh,ryh,rzh,bxh,byh,bzh
     &   ,mu,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mfdim,mfdim3,1)
      call mlw_gv(ro ,vx ,gx
     &  ,ee2,rx2,eeh,rxh
     &  ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,1)

      call rvtovv(vxh,vyh,vzh,rxh,ryh,rzh,roh,ix,jx,kx)
      call ee2eh(ehh,eeh,roh,vxh,vyh,vzh,bxh,byh,bzh,mu,ix,jx,kx)
      call eh2pr(prh,ehh,gm,ix,jx,kx)
      call ohmideal(exh,eyh,ezh,bxh,byh,bzh,vxh,vyh,vzh,ix,jx,kx)
      call mlw_mh_s(roh,prh,ehh,vxh,vyh,vzh,bxh,byh,bzh,exh,eyh,ezh
     &   ,xm,ym
     &   ,ro2,ee2,rx2,ry2,rz2,bx2,by2,bz2
     &   ,roh,eeh,rxh,ryh,rzh,bxh,byh,bzh
     &   ,mu,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mfdim,mfdim3,2)
      call mlw_gv(roh,vxh,gxm
     &  ,ee2,rx2,eeh,rxh
     &  ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,2)

      call rvtovv(vx,vy,vz,rx,ry,rz,ro,ix,jx,kx)
      qav=3.d0
      vvmin=1.d-4
      call avlap_8
     &    (ro2,ee2,rx2,ry2,rz2,bx2,by2,bz2
     &    ,ro ,ee ,rx ,ry ,rz ,bx ,by ,bz 
     &    ,vx,vy,vz,dt,qav,vvmin,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mfdim)

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
