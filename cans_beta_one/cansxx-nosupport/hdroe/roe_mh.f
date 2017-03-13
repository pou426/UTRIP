c======================================================================|
      subroutine roe_mh(rof,prf,vxf,vyf,vzf,bxf,byf,bzf
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

      dimension ron(ix,jx,kx),een(ix,jx,kx)
      dimension rxn(ix,jx,kx),ryn(ix,jx,kx),rzn(ix,jx,kx)
      dimension bxn(ix,jx,kx),byn(ix,jx,kx),bzn(ix,jx,kx)

      dimension fro(ix,jx,kx),fee(ix,jx,kx)
      dimension frx(ix,jx,kx),fry(ix,jx,kx),frz(ix,jx,kx)
      dimension fbx(ix,jx,kx),fby(ix,jx,kx),fbz(ix,jx,kx)

      dimension row(ix,jx,kx,2),prw(ix,jx,kx,2)
      dimension vxw(ix,jx,kx,2),vyw(ix,jx,kx,2),vzw(ix,jx,kx,2)
      dimension bxw(ix,jx,kx,2),byw(ix,jx,kx,2),bzw(ix,jx,kx,2)
c----------------------------------------------------------------------|

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
        call roeflux2_m(fro,fee,frx,fry,frz,fbx,fby,fbz,gm
     &               ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,ix,jx,kx)
      else if (mdir.eq.2) then
        call roeflux2_m(fro,fee,fry,frz,frx,fby,fbz,fbx,gm
     &               ,row,prw,vyw,vzw,vxw,byw,bzw,bxw,ix,jx,kx)
      else
        call roeflux2_m(fro,fee,frz,frx,fry,fbz,fbx,fby,gm
     &               ,row,prw,vzw,vxw,vyw,bzw,bxw,byw,ix,jx,kx)
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
