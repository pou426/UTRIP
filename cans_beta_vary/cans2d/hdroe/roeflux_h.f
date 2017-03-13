      subroutine roeflux_h(fro,fee,frx,fry,gm,row,prw,vxw,vyw,ix,jx)
      implicit double precision (a-h,o-z)

      dimension row(ix,jx,2),prw(ix,jx,2),vxw(ix,jx,2),vyw(ix,jx,2)
      dimension fro(ix,jx),fee(ix,jx),frx(ix,jx),fry(ix,jx)


      do j=1,jx-1
      do i=1,ix-1
         rhol=row(i,j,1)
         rhor=row(i,j,2)
         prl=prw(i,j,1)
         prr=prw(i,j,2)
         vxl=vxw(i,j,1)
         vxr=vxw(i,j,2)
         vyl=vyw(i,j,1)
         vyr=vyw(i,j,2)

c-----roe's variable
      sr0=sqrt(rhol)
      sr1=sqrt(rhor)
      sri=1.0d0/(sr0+sr1)
      rhobar=sr0*sr1
      vxbar=(sr0*vxl+sr1*vxr)*sri
      vybar=(sr0*vyl+sr1*vyr)*sri
      hl=0.5d0*(vxl**2+vyl**2)+gm*prl/((gm-1.0d0)*rhol)
      hr=0.5d0*(vxr**2+vyr**2)+gm*prr/((gm-1.0d0)*rhor)
      hbar=(sr0*hl+sr1*hr)*sri

c-----characteristic speed
      cs2=(gm-1.0d0)*(hbar-0.5d0*(vxbar**2+vybar**2))
      cs=sqrt(cs2)

c----- eigen value & entropy condition
      eeps=(vxr-vxl+abs(vxr-vxl))*2.5d-1 
      elpf=-max(abs(vxbar+cs),eeps)
      elmf=-max(abs(vxbar-cs),eeps)
      elze=-max(abs(vxbar),eeps)

c----- amplitude;w's
      drho=rhor-rhol
      du21=rhobar*(vxr-vxl)
      t2=(prr-prl)/cs

      wpf=0.5d0*(t2+du21)/cs
      wmf=0.5d0*(t2-du21)/cs
      wze=drho-(wpf+wmf)

c----- flux
      fluxlro=rhol*vxl
      fluxrro=rhor*vxr
      fluxlee=rhol*vxl*hl
      fluxree=rhor*vxr*hr
      fluxlrx=rhol*vxl**2+prl
      fluxrrx=rhor*vxr**2+prr
      fluxlry=rhol*vxl*vyl
      fluxrry=rhor*vxr*vyr

c----- components of the eigen vectors
      rpfro=1.
      rpfee=0.5d0*(vxbar**2+vybar**2)+cs*vxbar+cs2/(gm-1.0d0)
      rpfrx=vxbar+cs
      rpfry=vybar

      rmfro=1.
      rmfee=0.5d0*(vxbar**2+vybar**2)-cs*vxbar +cs2/(gm-1.0d0)
      rmfrx=vxbar-cs
      rmfry=vybar

      rzero=1.0d0
      rzeee=0.5d0*(vxbar**2+vybar**2)
      rzerx=vxbar
      rzery=vybar

c-----computation of f(i+1/2,j)
      fro(i,j)=0.5d0*(fluxlro+fluxrro
     1       +elpf*wpf*rpfro +elmf*wmf*rmfro
     2                                 +elze*wze*rzero)

      fee(i,j)=0.5d0*(fluxlee+fluxree
     1       +elpf*wpf*rpfee +elmf*wmf*rmfee 
     2                                 +elze*wze*rzeee)

      frx(i,j)=0.5d0*(fluxlrx+fluxrrx
     1       +elpf*wpf*rpfrx +elmf*wmf*rmfrx 
     2                                 +elze*wze*rzerx)

      fry(i,j)=0.5d0*(fluxlry+fluxrry
     1       +elpf*wpf*rpfry +elmf*wmf*rmfry 
     2                                 +elze*wze*rzery)

      enddo
      enddo

      return
      end
