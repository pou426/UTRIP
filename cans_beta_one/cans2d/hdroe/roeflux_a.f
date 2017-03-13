c======================================================================|
      subroutine roeflux_a(fro,row,vxw,ix,jx)
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
c    fro(ix,jx): [double] density flux
c
c INPUTS
c    row(ix,jx,2): [double] density at cell boundary
c    vxw(ix,jx,2): [double] velocity at cell boundary
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-7-1 T. Yokoyama based on N. Fukuda's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension row(ix,jx,2),vxw(ix,jx,2)
      dimension fro(ix,jx)
c----------------------------------------------------------------------|

      do j=1,jx-1
      do i=1,ix-1
        rhol=row(i,j,1)
        rhor=row(i,j,2)
        vxl=vxw(i,j,1)
        vxr=vxw(i,j,2)

c-----barred variables
        sr0=sqrt(rhol)
        sr1=sqrt(rhor)
        sri=1.0d0/(sr0+sr1)
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
        fro(i,j)=(fluxlro+fluxrro
     5            +elze*wze*rzero)*0.5d0

      enddo
      enddo

      return
      end
