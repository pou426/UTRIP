c======================================================================|
      subroutine roeflux_ht(fro,frx,fry,cs2,row,vxw,vyw,ix,jx)
c======================================================================|
c
c NAME  roeflux_ht
c
c PURPOSE
c    derive numerical flux by solving the linearized Riemann problem
c        * isothermal hydrodynamics
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
c    vxw(ix,2): [double] velocity at cell boundary
c    cs2: [double] square of sound speed
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on N. Fukuda's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension row(ix,jx,2),vxw(ix,jx,2),vyw(ix,jx,2)
      dimension fro(ix,jx),frx(ix,jx),fry(ix,jx)
c----------------------------------------------------------------------|

      do j=1,jx-1
      do i=1,ix-1
        rhol=row(i,j,1)
        rhor=row(i,j,2)
        vxl=vxw(i,j,1)
        vxr=vxw(i,j,2)
        vyl=vyw(i,j,1)
        vyr=vyw(i,j,2)
        prl=rhol*cs2
        prr=rhor*cs2

c-----barred variables
        sr0=sqrt(rhol)
        sr1=sqrt(rhor)
        sri=1.0d0/(sr0+sr1)
        rhobar=sr0*sr1
        vxbar=(sr0*vxl+sr1*vxr)*sri
        vybar=(sr0*vyl+sr1*vyr)*sri

c-----characteristic speed
        cs=sqrt(cs2)

c-----eigen values
        eeps=(vxr-vxl+abs(vxr-vxl))*2.5d-1 
        elpf=-max(abs(vxbar+cs),eeps)
        elmf=-max(abs(vxbar-cs),eeps)
        elze=-max(abs(vxbar),eeps)

c-----for singular points
        eps2=1.0d-6
        sgr2=cs2-eps2
        sp2=0.5d0+sign(0.5d0,sgr2)
        alphf=sp2*sqrt(cs2/(cs2+1.0d0-sp2))+1.0d0-sp2

c-----amplitude;w's
        drho=rhor-rhol
        du21=rhobar*(vxr-vxl)
        t2=drho*cs

        wpf=0.5d0*(t2+du21)/cs/alphf
        wmf=0.5d0*(t2-du21)/cs/alphf
        wze=drho-alphf*(wpf+wmf)

c-----components of the eigen vectors
        rpfro=alphf
        rpfrx=alphf*(vxbar+cs)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    need re-check by TY 2002-7-1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        rpfry=vybar

        rmfro=alphf
        rmfrx=alphf*(vxbar-cs)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    need re-check by TY 2002-7-1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        rmfry=vybar

        rzero=1.0d0
        rzerx=vxbar
        rzery=vybar

c----- flux_l, flux_r
        fluxlro=rhol*vxl
        fluxrro=rhor*vxr
        fluxlrx=rhol*vxl*vxl+prl
        fluxrrx=rhor*vxr*vxr+prr
        fluxlry=rhol*vxl*vyl
        fluxrry=rhor*vxr*vyr

c-----computation of f(i+1/2,j)
        fro(i,j)=(fluxlro+fluxrro
     1            +elpf*wpf*rpfro
     2            +elmf*wmf*rmfro
     5            +elze*wze*rzero)*0.5d0
        frx(i,j)=(fluxlrx+fluxrrx
     1            +elpf*wpf*rpfrx
     2            +elmf*wmf*rmfrx
     5            +elze*wze*rzerx)*0.5d0
        fry(i,j)=(fluxlry+fluxrry
     1            +elpf*wpf*rpfry
     2            +elmf*wmf*rmfry
     5            +elze*wze*rzery)*0.5d0

      enddo
      enddo

      return
      end
