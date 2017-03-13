c======================================================================|
      subroutine mlw_mh_s(rof,prf,ehf,vxf,vyf,vzf,bxf,byf,bzf
     &  ,exf,eyf,ezf,xf,yf
     &  ,ro2,ee2,rx2,ry2,rz2,bx2,by2,bz2
     &  ,roh,eeh,rxh,ryh,rzh,bxh,byh,bzh
     &  ,mu,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mfdim,mfdim3,mstage)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension mfdim(3),mfdim3(3,3)
      dimension dx(ix),dxm(ix),dxi(ix),dxim(ix),ux0(ix),ux1(ix)
      dimension dy(jx),dym(jx),dyi(jx),dyim(jx),uy0(jx),uy1(jx)
      dimension dz(kx),dzm(kx),dzi(kx),dzim(kx),uz0(kx),uz1(kx)

      dimension rof(ix,jx,kx),prf(ix,jx,kx),ehf(ix,jx,kx)
      dimension vxf(ix,jx,kx),vyf(ix,jx,kx),vzf(ix,jx,kx)
      dimension bxf(ix,jx,kx),byf(ix,jx,kx),bzf(ix,jx,kx)
      dimension exf(ix,jx,kx),eyf(ix,jx,kx),ezf(ix,jx,kx)
      dimension xf(ix),yf(jx)

      dimension roh(ix,jx,kx),eeh(ix,jx,kx)
      dimension rxh(ix,jx,kx),ryh(ix,jx,kx),rzh(ix,jx,kx)
      dimension bxh(ix,jx,kx),byh(ix,jx,kx),bzh(ix,jx,kx)

      dimension ro2(ix,jx,kx),ee2(ix,jx,kx)
      dimension rx2(ix,jx,kx),ry2(ix,jx,kx),rz2(ix,jx,kx)
      dimension bx2(ix,jx,kx),by2(ix,jx,kx),bz2(ix,jx,kx)

      dimension fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)
      dimension ss(ix,jx,kx)
      double precision mu
c----------------------------------------------------------------------|
c     ready
c----------------------------------------------------------------------|
      call dx2dxi(dx,dxi,ix)
      call dx2dxi(dxm,dxim,ix)
      call dx2ux(dx,dxm,ux0,ux1,ix)

      call dx2dxi(dy,dyi,jx)
      call dx2dxi(dym,dyim,jx)
      call dx2ux(dy,dym,uy0,uy1,jx)

      call dx2dxi(dz,dzi,kx)
      call dx2dxi(dzm,dzim,kx)
      call dx2ux(dz,dzm,uz0,uz1,kx)
c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      call getfro(fx,fy,fz,rof,vxf,vyf,vzf,ix,jx,kx,mfdim)
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss(i,j,k)=  fx(i,j,k)*2.d0+fy(i,j,k)/tan(yf(j))
         ss(i,j,k)= -ss(i,j,k)/xf(i)
         fy(i,j,k)= fy(i,j,k)/xf(i)
         fz(i,j,k)= fz(i,j,k)/xf(i)/sin(yf(j))
      enddo
      enddo
      enddo
      if (mstage.eq.1) then
        call mlwh(roh,ro2
     &      ,dt,fx,dxi,dxim,fy,dyi,dyim,fz,dzi,dzim,ix,jx,kx,mfdim)
        call mlwsh(roh,ro2,dt,ss,ix,jx,kx)
      else
        call mlwsf(ro2,dt,ss,ux0,ux1,uy0,uy1,uz0,uz1,ix,jx,kx)
        call mlwf(ro2
     &       ,dt,fx,dxi,ux0,ux1,fy,dyi,uy0,uy1,fz,dzi,uz0,uz1
     &       ,ix,jx,kx,mfdim)
      endif

c---  energy ---
      call getfee(fx,fy,fz
     &      ,rof,prf,ehf,vxf,vyf,vzf,bxf,byf,bzf,exf,eyf,ezf
     &      ,mu,ix,jx,kx,mfdim)
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss(i,j,k)=  fx(i,j,k)*2.d0+fy(i,j,k)/tan(yf(j))
         ss(i,j,k)= -ss(i,j,k)/xf(i)
         fy(i,j,k)= fy(i,j,k)/xf(i)
         fz(i,j,k)= fz(i,j,k)/xf(i)/sin(yf(j))
      enddo
      enddo
      enddo
      if (mstage.eq.1) then
        call mlwh(eeh,ee2 
     &      ,dt,fx,dxi,dxim,fy,dyi,dyim,fz,dzi,dzim,ix,jx,kx,mfdim)
        call mlwsh(eeh,ee2,dt,ss,ix,jx,kx)
      else
        call mlwsf(ee2,dt,ss,ux0,ux1,uy0,uy1,uz0,uz1,ix,jx,kx)
        call mlwf(ee2
     &       ,dt,fx,dxi,ux0,ux1,fy,dyi,uy0,uy1,fz,dzi,uz0,uz1
     &       ,ix,jx,kx,mfdim)
      endif

c---  x-momentum ---
      call getfrx(fx,fy,fz,rof,prf,vxf,vyf,vzf,bxf,byf,bzf,mu
     &      ,ix,jx,kx,mfdim3(:,1))
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss(i,j,k)= 
     &         +(vxf(i,j,k)*vxf(i,j,k)*rof(i,j,k)
     &          -bxf(i,j,k)*bxf(i,j,k)/mu)*2.d0
     &         +(vxf(i,j,k)*vyf(i,j,k)*rof(i,j,k)
     &          -bxf(i,j,k)*byf(i,j,k)/mu)/tan(yf(j))
     &         -(vyf(i,j,k)*vyf(i,j,k)*rof(i,j,k)
     &          -byf(i,j,k)*byf(i,j,k)/mu)
     &         -(vzf(i,j,k)*vzf(i,j,k)*rof(i,j,k)
     &          -bzf(i,j,k)*bzf(i,j,k)/mu)
         ss(i,j,k)= -ss(i,j,k)/xf(i)
         fy(i,j,k)= fy(i,j,k)/xf(i)
         fz(i,j,k)= fz(i,j,k)/xf(i)/sin(yf(j))
      enddo
      enddo
      enddo
      if (mstage.eq.1) then
        call mlwh(rxh,rx2
     &    ,dt,fx,dxi,dxim,fy,dyi,dyim,fz,dzi,dzim,ix,jx,kx,mfdim)
        call mlwsh(rxh,rx2,dt,ss,ix,jx,kx)
      else
        call mlwsf(rx2,dt,ss,ux0,ux1,uy0,uy1,uz0,uz1,ix,jx,kx)
        call mlwf(rx2
     &       ,dt,fx,dxi,ux0,ux1,fy,dyi,uy0,uy1,fz,dzi,uz0,uz1
     &       ,ix,jx,kx,mfdim)
      endif

c---  y-momentum ---
      call getfrx(fy,fz,fx,rof,prf,vyf,vzf,vxf,byf,bzf,bxf,mu
     &      ,ix,jx,kx,mfdim3(:,2))
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss(i,j,k)= 
     &         +(vyf(i,j,k)*vxf(i,j,k)*rof(i,j,k)
     &          -byf(i,j,k)*bxf(i,j,k)/mu)
     &         +(vyf(i,j,k)*vyf(i,j,k)*rof(i,j,k)
     &          -byf(i,j,k)*byf(i,j,k)/mu)/tan(yf(j))
     &         -(vzf(i,j,k)*vzf(i,j,k)*rof(i,j,k)
     &          -bzf(i,j,k)*bzf(i,j,k)/mu)/tan(yf(j))
     &         +(vxf(i,j,k)*vzf(i,j,k)*rof(i,j,k)
     &          -bxf(i,j,k)*bzf(i,j,k)/mu)*2.d0
         ss(i,j,k)= -ss(i,j,k)/xf(i)
         fy(i,j,k)= fy(i,j,k)/xf(i)
         fz(i,j,k)= fz(i,j,k)/xf(i)/sin(yf(j))
      enddo
      enddo
      enddo
      if (mstage.eq.1) then
        call mlwh(ryh,ry2
     &    ,dt,fx,dxi,dxim,fy,dyi,dyim,fz,dzi,dzim,ix,jx,kx,mfdim)
        call mlwsh(ryh,ry2,dt,ss,ix,jx,kx)
      else
        call mlwsf(ry2,dt,ss,ux0,ux1,uy0,uy1,uz0,uz1,ix,jx,kx)
        call mlwf(ry2
     &       ,dt,fx,dxi,ux0,ux1,fy,dyi,uy0,uy1,fz,dzi,uz0,uz1
     &       ,ix,jx,kx,mfdim)
      endif
c      do j=30,40
c      write(6,*) ryh(5,j,2)
c      enddo
c      stop

c---  z-momentum ---
      call getfrx(fz,fx,fy,rof,prf,vzf,vxf,vyf,bzf,bxf,byf,mu
     &      ,ix,jx,kx,mfdim3(:,3))
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss(i,j,k)= 
     &         +(vzf(i,j,k)*vxf(i,j,k)*rof(i,j,k)
     &          -bzf(i,j,k)*bxf(i,j,k)/mu)*3.d0
     &         +(vzf(i,j,k)*vyf(i,j,k)*rof(i,j,k)
     &          -bzf(i,j,k)*byf(i,j,k)/mu)*2.d0/tan(yf(j))
         ss(i,j,k)= -ss(i,j,k)/xf(i)
         fy(i,j,k)= fy(i,j,k)/xf(i)
         fz(i,j,k)= fz(i,j,k)/xf(i)/sin(yf(j))
      enddo
      enddo
      enddo
      if (mstage.eq.1) then
        call mlwh(rzh,rz2
     &    ,dt,fx,dxi,dxim,fy,dyi,dyim,fz,dzi,dzim,ix,jx,kx,mfdim)
        call mlwsh(rzh,rz2,dt,ss,ix,jx,kx)
      else
        call mlwsf(rz2,dt,ss,ux0,ux1,uy0,uy1,uz0,uz1,ix,jx,kx)
        call mlwf(rz2,dt
     &       ,fx,dxi,ux0,ux1,fy,dyi,uy0,uy1,fz,dzi,uz0,uz1
     &       ,ix,jx,kx,mfdim)
      endif

c---  x-magnetic ---
      call getfbx(fx,fy,fz,eyf,ezf,ix,jx,kx,mfdim3(:,1))
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss(i,j,k)= ezf(i,j,k)/tan(yf(j))
         ss(i,j,k)= -ss(i,j,k)/xf(i)
         fy(i,j,k)= fy(i,j,k)/xf(i)
         fz(i,j,k)= fz(i,j,k)/xf(i)/sin(yf(j))
      enddo
      enddo
      enddo
      if (mstage.eq.1) then
        call mlwh(bxh,bx2
     &    ,dt,fx,dxi,dxim,fy,dyi,dyim,fz,dzi,dzim,ix,jx,kx,mfdim)
        call mlwsh(bxh,bx2,dt,ss,ix,jx,kx)
      else
        call mlwsf(bx2,dt,ss,ux0,ux1,uy0,uy1,uz0,uz1,ix,jx,kx)
        call mlwf(bx2,dt
     &       ,fx,dxi,ux0,ux1,fy,dyi,uy0,uy1,fz,dzi,uz0,uz1
     &       ,ix,jx,kx,mfdim)
      endif

c---  y-magnetic ---
      call getfbx(fy,fz,fx,ezf,exf,ix,jx,kx,mfdim3(:,2))
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss(i,j,k)=-ezf(i,j,k)
         ss(i,j,k)= -ss(i,j,k)/xf(i)
         fy(i,j,k)= fy(i,j,k)/xf(i)
         fz(i,j,k)= fz(i,j,k)/xf(i)/sin(yf(j))
      enddo
      enddo
      enddo
      if (mstage.eq.1) then
        call mlwh(byh,by2
     &    ,dt,fx,dxi,dxim,fy,dyi,dyim,fz,dzi,dzim,ix,jx,kx,mfdim)
        call mlwsh(byh,by2,dt,ss,ix,jx,kx)
      else
        call mlwsf(by2,dt,ss,ux0,ux1,uy0,uy1,uz0,uz1,ix,jx,kx)
        call mlwf(by2
     &       ,dt,fx,dxi,ux0,ux1,fy,dyi,uy0,uy1,fz,dzi,uz0,uz1
     &       ,ix,jx,kx,mfdim)
      endif

c---  z-magnetic ---
      call getfbx(fz,fx,fy,exf,eyf,ix,jx,kx,mfdim3(:,3))
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss(i,j,k)= eyf(i,j,k)
         ss(i,j,k)= -ss(i,j,k)/xf(i)
         fy(i,j,k)= fy(i,j,k)/xf(i)
         fz(i,j,k)= fz(i,j,k)/xf(i)/sin(yf(j))
      enddo
      enddo
      enddo
      if (mstage.eq.1) then
        call mlwh(bzh,bz2
     &    ,dt,fx,dxi,dxim,fy,dyi,dyim,fz,dzi,dzim,ix,jx,kx,mfdim)
        call mlwsh(bzh,bz2,dt,ss,ix,jx,kx)
      else
        call mlwsf(bz2,dt,ss,ux0,ux1,uy0,uy1,uz0,uz1,ix,jx,kx)
        call mlwf(bz2
     &       ,dt,fx,dxi,ux0,ux1,fy,dyi,uy0,uy1,fz,dzi,uz0,uz1
     &       ,ix,jx,kx,mfdim)
      endif


      return
      end
