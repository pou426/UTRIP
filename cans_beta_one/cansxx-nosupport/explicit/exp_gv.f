c======================================================================|
      subroutine exp_gv(rxn,een,rof,vxf,gxf,dt,ix,jx,kx)
c======================================================================|
      implicit double precision (a-h,o-z)

      dimension rof(ix,jx,kx),vxf(ix,jx,kx)
      dimension een(ix,jx,kx),rxn(ix,jx,kx)

      dimension gxf(ix,jx,kx)

c----------------------------------------------------------------------|
c     gravity
c----------------------------------------------------------------------|
c---  x-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss= rof(i,j,k)*gxf(i,j,k)
         rxn(i,j,k)=rxn(i,j,k)+dt*ss
      enddo
      enddo
      enddo

c---  x-comp. energy ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss= rof(i,j,k)*vxf(i,j,k)*gxf(i,j,k)
         een(i,j,k)=een(i,j,k)+dt*ss
      enddo
      enddo
      enddo


      return
      end
