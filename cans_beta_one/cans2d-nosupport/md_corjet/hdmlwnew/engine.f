      subroutine engine(ro,pr,vx,vy,bx,by,az,dt,gm
     &               ,gx,gxm,gy,gym,et,etm,dx,dxm,dy,dym,ix,jx)

      implicit double precision (a-h,o-z)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
     &         ,bx(ix,jx),by(ix,jx),az(ix,jx)
      dimension gx(ix,jx),gxm(ix,jx)
      dimension gy(ix,jx),gym(ix,jx)
      dimension et(ix,jx),etm(ix,jx)

      dimension ee(ix,jx),rx(ix,jx),ry(ix,jx)
      dimension roh(ix,jx),eeh(ix,jx),rxh(ix,jx),ryh(ix,jx)
     &         ,bxh(ix,jx),byh(ix,jx),azh(ix,jx)
      dimension ro2(ix,jx),ee2(ix,jx),rx2(ix,jx),ry2(ix,jx)
     &         ,bx2(ix,jx),by2(ix,jx),az2(ix,jx)
c----------------------------------------------------------------------|


      call rv2vv(ro,rx,ry,vx,vy,ix,jx,-1)
      call ee2pr(pr,ee,ro,rx,ry,bx,by,gm,ix,jx,-1)
      call mlwreset(ro2,roh,ro,ix,jx)
      call mlwreset(ee2,eeh,ee,ix,jx)
      call mlwreset(rx2,rxh,rx,ix,jx)
      call mlwreset(ry2,ryh,ry,ix,jx)
      call mlwreset(bx2,bxh,bx,ix,jx)
      call mlwreset(by2,byh,by,ix,jx)
      call mlwreset(az2,azh,az,ix,jx)

      call mlw_mh(ro,ee,rx,ry,bx,by
     &  ,ro2,ee2,rx2,ry2,bx2,by2,az2
     &  ,roh,eeh,rxh,ryh,bxh,byh,azh
     &  ,gm,dt,dx,dxm,dy,dym,ix,jx,1)
      call mlw_e(bx,by,et,bx2,by2,az2,bxh,byh,azh
     &  ,dt,dx,dxm,dy,dym,ix,jx,1)
c     call mlw_g(ro,ee,rx,gx,ee2,ry2,eeh,rxh
c    &  ,dt,dx,dxm,dy,dym,ix,jx,1)
      call mlw_g(ro,ee,ry,gy,ee2,ry2,eeh,ryh
     &  ,dt,dx,dxm,dy,dym,ix,jx,1)

      call mlw_mh(roh,eeh,rxh,ryh,bxh,byh
     &  ,ro2,ee2,rx2,ry2,bx2,by2,az2
     &  ,roh,eeh,rxh,ryh,bxh,byh,azh
     &  ,gm,dt,dx,dxm,dy,dym,ix,jx,2)
      call mlw_e(bxh,byh,etm,bx2,by2,az2,bxh,byh,azh
     &  ,dt,dx,dxm,dy,dym,ix,jx,2)
c     call mlw_g(roh,eeh,rxh,gxm,ee2,ry2,eeh,rxh
c    &  ,dt,dx,dxm,dy,dym,ix,jx,1)
      call mlw_g(roh,eeh,ryh,gym,ee2,ry2,eeh,ryh
     &  ,dt,dx,dxm,dy,dym,ix,jx,1)

      call rv2vv(ro,rx,ry,vx,vy,ix,jx,+1)
      qav=3.d0
      vvmin=1.d-4
      call avlap_6(ro2,ee2,rx2,ry2,bx2,by2,ro,ee,rx,ry,bx,by
     &    ,vx,vy,dt,qav,vvmin,dx,dxm,dy,dym,ix,jx)

      call daupdate(ro,ro2,ix,jx)
      call daupdate(ee,ee2,ix,jx)
      call daupdate(rx,rx2,ix,jx)
      call daupdate(ry,ry2,ix,jx)
      call daupdate(bx,bx2,ix,jx)
      call daupdate(by,by2,ix,jx)
      call daupdate(az,az2,ix,jx)

      call rv2vv(ro,rx,ry,vx,vy,ix,jx,+1)
      call ee2pr(pr,ee,ro,rx,ry,bx,by,gm,ix,jx,+1)

      return
      end
