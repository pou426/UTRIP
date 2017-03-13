c======================================================================|
      subroutine mlw_gv(rof,vxf,gxf
     &  ,ee2,rx2,eeh,rxh
     &  ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mfdim,mstage)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension mfdim(3)
      dimension dx(ix),dxm(ix),ux0(ix),ux1(ix)
      dimension dy(jx),dym(jx),uy0(jx),uy1(jx)
      dimension dz(kx),dzm(kx),uz0(kx),uz1(kx)

      dimension rof(ix,jx,kx),vxf(ix,jx,kx)
      dimension ee2(ix,jx,kx),rx2(ix,jx,kx)
      dimension eeh(ix,jx,kx),rxh(ix,jx,kx)

      dimension gxf(ix,jx,kx)

      dimension ss(ix,jx,kx)
c----------------------------------------------------------------------|
c     ready
c----------------------------------------------------------------------|
      call dx2ux(dx,dxm,ux0,ux1,ix)
      call dx2ux(dy,dym,uy0,uy1,jx)
      call dx2ux(dz,dzm,uz0,uz1,kx)
c----------------------------------------------------------------------|
c     gravity
c----------------------------------------------------------------------|
c---  x-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss(i,j,k)= rof(i,j,k)*gxf(i,j,k)
      enddo
      enddo
      enddo
      if (mstage.eq.1) then
        call mlwsh(rxh,rx2,dt,ss,ix,jx,kx,mfdim)
      else
        call mlwsf(rx2,dt,ss,ux0,ux1,uy0,uy1,uz0,uz1,ix,jx,kx,mfdim)
      endif

c---  x-comp. energy ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss(i,j,k)= rof(i,j,k)*vxf(i,j,k)*gxf(i,j,k)
      enddo
      enddo
      enddo
      if (mstage.eq.1) then
        call mlwsh(eeh,ee2,dt,ss,ix,jx,kx,mfdim)
      else
        call mlwsf(ee2,dt,ss,ux0,ux1,uy0,uy1,uz0,uz1,ix,jx,kx,mfdim)
      endif


      return
      end
