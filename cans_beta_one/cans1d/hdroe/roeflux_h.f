c======================================================================|
      subroutine roeflux_h(fro,fee,frx,gm,row,prw,vxw,ix)
c======================================================================|
c
c NAME  roeflux_h
c
c PURPOSE
c    derive numerical flux by solving the linearized Riemann problem
c        * hydrodynamics
c
c INPUTS & OUTPUTS
c    None
c
c OUTPUTS
c    fro(ix): [double] density flux
c    fee(ix): [double] total-energy flux
c    frx(ix): [double] momentum flux
c
c INPUTS
c    row(ix,2): [double] density at cell boundary
c    prw(ix,2): [double] pressure at cell boundary
c    vxw(ix,2): [double] velocity at cell boundary
c    gm: [double] polytropic index gamma
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on N. Fukuda's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension row(ix,2),prw(ix,2),vxw(ix,2)
      dimension fro(ix),fee(ix),frx(ix)
c----------------------------------------------------------------------|


      do i=1,ix-1
         rhol=row(i,1)
         vxl=vxw(i,1)
         prl=prw(i,1)
         rhor=row(i,2)
         vxr=vxw(i,2)
         prr=prw(i,2)

c-----roe's variable
      sr0=sqrt(rhol)
      sr1=sqrt(rhor)
      sri=1.0d0/(sr0+sr1)
      rhobar=sr0*sr1
      vxbar=(sr0*vxl+sr1*vxr)*sri
      hl=0.5d0*(vxl**2)+gm*prl/((gm-1.0d0)*rhol)
      hr=0.5d0*(vxr**2)+gm*prr/((gm-1.0d0)*rhor)
      hbar=(sr0*hl+sr1*hr)*sri

c-----characteristic speed
      cs2=(gm-1.0d0)*(hbar-0.5d0*(vxbar**2))
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

c----- components of the eigen vectors
      rpfro=1.
      rpfrx=vxbar+cs
      rpfee=0.5d0*vxbar**2+cs*vxbar+cs2/(gm-1.0d0)

      rmfro=1.
      rmfrx=vxbar-cs
      rmfee=0.5d0*vxbar**2-cs*vxbar +cs2/(gm-1.0d0)

      rzero=1.0d0
      rzerx=vxbar
      rzeee=0.5d0*vxbar**2

c----- flux_l
      fluxlro=rhol*vxl
      fluxlrx=rhol*vxl**2+prl
      fluxlee=rhol*vxl*hl
c----- flux_r
      fluxrro=rhor*vxr
      fluxrrx=rhor*vxr**2+prr
      fluxree=rhor*vxr*hr

c-----computation of f(i+1/2,j)
      fro(i)=0.5d0*(fluxlro+fluxrro
     1       +elpf*wpf*rpfro +elmf*wmf*rmfro
     2                                 +elze*wze*rzero)

      frx(i)=0.5d0*(fluxlrx+fluxrrx
     1       +elpf*wpf*rpfrx +elmf*wmf*rmfrx 
     2                                 +elze*wze*rzerx)

      fee(i)=0.5d0*(fluxlee+fluxree
     1       +elpf*wpf*rpfee +elmf*wmf*rmfee 
     2                                 +elze*wze*rzeee)

      enddo

      return
      end
