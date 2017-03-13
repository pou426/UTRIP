c======================================================================|
      subroutine exp_et
     &    (bxn,byn,bzn,een
     &    ,bxf,byf,bzf
     &    ,et,mu,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mfdim,mfdim3)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension mfdim(3),mfdim3(3,3)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension dz(kx),dzm(kx)

      dimension bxn(ix,jx,kx),byn(ix,jx,kx),bzn(ix,jx,kx),een(ix,jx,kx)
      dimension bxf(ix,jx,kx),byf(ix,jx,kx),bzf(ix,jx,kx)

      dimension cxf(ix,jx,kx),cyf(ix,jx,kx),czf(ix,jx,kx)
      dimension exf(ix,jx,kx),eyf(ix,jx,kx),ezf(ix,jx,kx)

      dimension et(ix,jx,kx)

      dimension fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)

      double precision mu

c----------------------------------------------------------------------|

      call bbtocc(cxf,cyf,czf,bxf,byf,bzf,dx,dy,dz,ix,jx,kx,mfdim)

      do k=1,kx
      do j=1,jx
      do i=1,ix
         exf(i,j,k) = et(i,j,k)*cxf(i,j,k)
         eyf(i,j,k) = et(i,j,k)*cyf(i,j,k)
         ezf(i,j,k) = et(i,j,k)*czf(i,j,k)
      enddo
      enddo
      enddo

c---  x-magnetic ---
      call getfbx(fx,fy,fz,eyf,ezf,ix,jx,kx,mfdim3(:,1))
      if (mfdim(1).eq.1) call adddfdx(bxn,fx,dt,dx,dy,dz,ix,jx,kx,1)
      if (mfdim(2).eq.1) call adddfdx(bxn,fy,dt,dx,dy,dz,ix,jx,kx,2)
      if (mfdim(3).eq.1) call adddfdx(bxn,fz,dt,dx,dy,dz,ix,jx,kx,3)

c---  y-magnetic ---
      call getfbx(fy,fz,fx,ezf,exf,ix,jx,kx,mfdim3(:,2))
      if (mfdim(1).eq.1) call adddfdx(byn,fx,dt,dx,dy,dz,ix,jx,kx,1)
      if (mfdim(2).eq.1) call adddfdx(byn,fy,dt,dx,dy,dz,ix,jx,kx,2)
      if (mfdim(3).eq.1) call adddfdx(byn,fz,dt,dx,dy,dz,ix,jx,kx,3)

c---  z-magnetic ---
      call getfbx(fz,fx,fy,exf,eyf,ix,jx,kx,mfdim3(:,3))
      if (mfdim(1).eq.1) call adddfdx(bzn,fx,dt,dx,dy,dz,ix,jx,kx,1)
      if (mfdim(2).eq.1) call adddfdx(bzn,fy,dt,dx,dy,dz,ix,jx,kx,2)
      if (mfdim(3).eq.1) call adddfdx(bzn,fz,dt,dx,dy,dz,ix,jx,kx,3)

c---  energy ---
      call getfeb(fx,fy,fz,bxf,byf,bzf,exf,eyf,ezf
     &      ,mu,ix,jx,kx,mfdim)
      if (mfdim(1).eq.1) call adddfdx(een,fx,dt,dx,dy,dz,ix,jx,kx,1)
      if (mfdim(2).eq.1) call adddfdx(een,fy,dt,dx,dy,dz,ix,jx,kx,2)
      if (mfdim(3).eq.1) call adddfdx(een,fz,dt,dx,dy,dz,ix,jx,kx,3)

      return
      end
