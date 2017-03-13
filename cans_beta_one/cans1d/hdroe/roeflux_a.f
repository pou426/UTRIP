c======================================================================|
      subroutine roeflux_a(fro,row,vxw,ix)
c======================================================================|
c
c NAME  roeflux_a
c
c PURPOSE
c    derive numerical flux by solving the linearized Riemann problem
c        * simple advection
c
c INPUTS & OUTPUTS
c    None
c
c OUTPUTS
c    fro(ix): [double] density flux
c
c INPUTS
c    row(ix,2): [double] density at cell boundary
c    vxw(ix,2): [double] velocity at cell boundary
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on N. Fukuda's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension row(ix,2),vxw(ix,2)
      dimension fro(ix)
c----------------------------------------------------------------------|

      do i=1,ix-1
        rhol=row(i,1)
        vxl=vxw(i,1)
        rhor=row(i,2)
        vxr=vxw(i,2)

c-----barred variables
        sr0=sqrt(rhol)
        sr1=sqrt(rhor)
        sri=1.0d0/(sr0+sr1)
        rhobar=sr0*sr1
        vxbar=(sr0*vxl+sr1*vxr)*sri

c-----eigen values
        eeps=(vxr-vxl+abs(vxr-vxl))*2.5d-1 
        elze=-max(abs(vxbar),eeps)

c-----amplitude;w's
        wze=rhor-rhol

c-----components of the eigen vectors

        rzero=1.0d0

c----- flux_l, flux_r
        fluxlro=rhol*vxl
        fluxrro=rhor*vxr

c-----computation of f(i+1/2,j)
        fro(i)=(fluxlro+fluxrro
     5            +elze*wze*rzero)*0.5d0

      enddo

      return
      end
