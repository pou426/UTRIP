c======================================================================|
      subroutine mlw_et(bxf,byf,bzf,etf
     &  ,bx2,by2,bz2,bxh,byh,bzh
     &  ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mfdim,mfdim3,mstage)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension mfdim(3),mfdim3(3,3)
      dimension dx(ix),dxm(ix),dxi(ix),dxim(ix),ux0(ix),ux1(ix)
      dimension dy(jx),dym(jx),dyi(jx),dyim(jx),uy0(jx),uy1(jx)
      dimension dz(kx),dzm(kx),dzi(kx),dzim(kx),uz0(kx),uz1(kx)

      dimension bxf(ix,jx,kx),byf(ix,jx,kx),bzf(ix,jx,kx)
      dimension exf(ix,jx,kx),eyf(ix,jx,kx),ezf(ix,jx,kx)
      dimension cxf(ix,jx,kx),cyf(ix,jx,kx),czf(ix,jx,kx)

      dimension bxh(ix,jx,kx),byh(ix,jx,kx),bzh(ix,jx,kx)
      dimension bx2(ix,jx,kx),by2(ix,jx,kx),bz2(ix,jx,kx)

      dimension fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)
      dimension etf(ix,jx,kx)
c----------------------------------------------------------------------|

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

      call bbtocx(cxf,byf,bzf,dy,dz,ix,jx,kx)
      call bbtocy(cyf,bzf,bxf,dz,dx,ix,jx,kx)
      call bbtocz(czf,bxf,byf,dx,dy,ix,jx,kx)

      do k=1,kx
      do j=1,jx
      do i=1,ix
         exf(i,j,k) =  etf(i,j,k)*cxf(i,j,k)
         eyf(i,j,k) =  etf(i,j,k)*cyf(i,j,k)
         ezf(i,j,k) =  etf(i,j,k)*czf(i,j,k)
      enddo
      enddo
      enddo

c---  x-magnetic ---
      call getfbx(fx,fy,fz,eyf,ezf,ix,jx,kx,mfdim3(:,1))
      if (mstage.eq.1) then
        call mlwh(bxh,bx2
     &    ,dt,fx,dxi,dxim,fy,dyi,dyim,fz,dzi,dzim,ix,jx,kx,mfdim)
      else
        call mlwf(bx2,dt
     &       ,fx,dxi,ux0,ux1,fy,dyi,uy0,uy1,fz,dzi,uz0,uz1
     &       ,ix,jx,kx,mfdim)
      endif

c---  y-magnetic ---
      call getfbx(fy,fz,fx,ezf,exf,ix,jx,kx,mfdim3(:,2))
      if (mstage.eq.1) then
        call mlwh(byh,by2
     &    ,dt,fx,dxi,dxim,fy,dyi,dyim,fz,dzi,dzim,ix,jx,kx,mfdim)
      else
        call mlwf(by2
     &       ,dt,fx,dxi,ux0,ux1,fy,dyi,uy0,uy1,fz,dzi,uz0,uz1
     &       ,ix,jx,kx,mfdim)
      endif

c---  z-magnetic ---
      call getfbx(fz,fx,fy,exf,eyf,ix,jx,kx,mfdim3(:,3))
      if (mstage.eq.1) then
        call mlwh(bzh,bz2
     &    ,dt,fx,dxi,dxim,fy,dyi,dyim,fz,dzi,dzim,ix,jx,kx,mfdim)
      else
        call mlwf(bz2
     &       ,dt,fx,dxi,ux0,ux1,fy,dyi,uy0,uy1,fz,dzi,uz0,uz1
     &       ,ix,jx,kx,mfdim)
      endif


      return
      end
