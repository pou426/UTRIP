c======================================================================|
      subroutine hlld_mh(rof,prf,vxf,vyf,vzf,bxf,byf,bzf
     &   ,ron,een,rxn,ryn,rzn,bxn,byn,bzn
     &   ,gm,dt,dx,dy,dz,ix,jx,kx,mfdim,mstage)
c======================================================================|
c     numerical solver of mhd equations by roe method with muscl
c     for ideal 1d simulation (2nd order)
c     version 1.1 (2001/08/24 naoya fukuda)
c     version 1.2 (2005/02/08 Yuji SATO)
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)

      dimension mfdim(3)
      dimension dx(ix),dy(jx),dz(kx)

      dimension rof(ix,jx,kx),prf(ix,jx,kx)
      dimension vxf(ix,jx,kx),vyf(ix,jx,kx),vzf(ix,jx,kx)
      dimension bxf(ix,jx,kx),byf(ix,jx,kx),bzf(ix,jx,kx)

      dimension bb(ix,jx,kx)
      dimension ron(ix,jx,kx),een(ix,jx,kx)
      dimension rxn(ix,jx,kx),ryn(ix,jx,kx),rzn(ix,jx,kx)
      dimension bxn(ix,jx,kx),byn(ix,jx,kx),bzn(ix,jx,kx)

      dimension fro(ix,jx,kx),fee(ix,jx,kx)
      dimension frx(ix,jx,kx),fry(ix,jx,kx),frz(ix,jx,kx)
      dimension fbx(ix,jx,kx),fby(ix,jx,kx),fbz(ix,jx,kx)

      dimension row(ix,jx,kx,2),prw(ix,jx,kx,2)
      dimension vxw(ix,jx,kx,2),vyw(ix,jx,kx,2),vzw(ix,jx,kx,2)
      dimension bxw(ix,jx,kx,2),byw(ix,jx,kx,2),bzw(ix,jx,kx,2)
      double precision mu
c----------------------------------------------------------------------|
      floor=1.d-10

      do mdir=1,3

      if (mfdim(mdir).eq.1) then

      if (mstage.eq.1) then
      call inp_const(rof,row,ix,jx,kx,mdir)
      call inp_const(prf,prw,ix,jx,kx,mdir)
      call inp_const(vxf,vxw,ix,jx,kx,mdir)
      call inp_const(vyf,vyw,ix,jx,kx,mdir)
      call inp_const(vzf,vzw,ix,jx,kx,mdir)
      call inp_const(bxf,bxw,ix,jx,kx,mdir)
      call inp_const(byf,byw,ix,jx,kx,mdir)
      call inp_const(bzf,bzw,ix,jx,kx,mdir)
      else
      call inp_tvdminmod(rof,row,ix,jx,kx,mdir)
      call inp_tvdminmod(prf,prw,ix,jx,kx,mdir)
      call inp_tvdminmod(vxf,vxw,ix,jx,kx,mdir)
      call inp_tvdminmod(vyf,vyw,ix,jx,kx,mdir)
      call inp_tvdminmod(vzf,vzw,ix,jx,kx,mdir)
      call inp_tvdminmod(bxf,bxw,ix,jx,kx,mdir)
      call inp_tvdminmod(byf,byw,ix,jx,kx,mdir)
      call inp_tvdminmod(bzf,bzw,ix,jx,kx,mdir)
      endif

      if (mdir.eq.1) then
        do k=1,kx
        do j=1,jx
        do i=1,ix-1
          bb(i,j,k)=(bxf(i,j,k)+bxf(i+1,j,k))/2.d0
        enddo
        enddo
        enddo
        floor=1.d-10
        call hlld_flux(row,prw,vxw,vyw,vzw,bb,byw,bzw,gm,ix,jx,kx 
     &               ,floor,fro,fee,frx,fry,frz,fby,fbz)
        do k=1,kx
        do j=1,jx
        do i=1,ix-1
          fbx(i,j,k)=0.d0
        enddo
        enddo
        enddo
      else if (mdir.eq.2) then
        do k=1,kx
        do j=1,jx-1
        do i=1,ix
          bb(i,j,k)=(byf(i,j,k)+byf(i,j+1,k))/2.d0
        enddo
        enddo
        enddo
        call hlld_flux(row,prw,vyw,vzw,vxw,bb,bzw,bxw,gm,ix,jx,kx 
     &               ,floor,fro,fee,fry,frz,frx,fbz,fbx)
        do k=1,kx
        do j=1,jx
        do i=1,ix-1
          fby(i,j,k)=0.d0
        enddo
        enddo
        enddo
      else
        do k=1,kx-1
        do j=1,jx
        do i=1,ix
          bb(i,j,k)=(bzf(i,j,k)+bzf(i,j,k+1))/2.d0
        enddo
        enddo
        enddo
        call hlld_flux(row,prw,vzw,vxw,vyw,bb,bxw,byw,gm,ix,jx,kx 
     &               ,floor,fro,fee,frz,frx,fry,fbx,fby)
        do k=1,kx
        do j=1,jx
        do i=1,ix-1
          fbz(i,j,k)=0.d0
        enddo
        enddo
        enddo
      endif

      call adddfdx(ron,fro,dt,dx,dy,dz,ix,jx,kx,mdir)
      call adddfdx(een,fee,dt,dx,dy,dz,ix,jx,kx,mdir)
      call adddfdx(rxn,frx,dt,dx,dy,dz,ix,jx,kx,mdir)
      call adddfdx(ryn,fry,dt,dx,dy,dz,ix,jx,kx,mdir)
      call adddfdx(rzn,frz,dt,dx,dy,dz,ix,jx,kx,mdir)
      call adddfdx(bxn,fbx,dt,dx,dy,dz,ix,jx,kx,mdir)
      call adddfdx(byn,fby,dt,dx,dy,dz,ix,jx,kx,mdir)
      call adddfdx(bzn,fbz,dt,dx,dy,dz,ix,jx,kx,mdir)
      endif

      enddo

      return
      end
