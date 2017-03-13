c======================================================================|
      subroutine ctransptz(bzh,vznx,vzny,bznx,bzny,vxm,vym,bxm,bym
     &      ,dt,dx,dy,ix,jx)
c======================================================================|
c
c NAME  moc
c
c PURPOSE
c    Constrained Transport (CT) method for induction equation
c    satisfying div B = 0.
c
c INPUTS & OUTPUTS
c    byh(ix): [double] magnetic field
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    vym(ix): [double] velocity
c    bym(ix): [double] magnetic field
c    vxm(ix) : [double] velocity along the x-cordinate
c    bxm(ix) : [double] magnetic field
c    dx(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dx(ix),dy(jx)

      dimension vznx(ix,jx),vzny(ix,jx),bznx(ix,jx),bzny(ix,jx)
      dimension vxm(ix,jx),vym(ix,jx),bxm(ix,jx),bym(ix,jx)
      dimension eynx(ix,jx),exny(ix,jx)
      dimension bzh(ix,jx)

c----------------------------------------------------------------------|


      do j=1,jx
      do i=1,ix-1
        eynx(i,j)=vxm(i,j)*bznx(i,j)-vznx(i,j)*bxm(i,j)
      enddo
      enddo
      do j=1,jx
      do i=1,ix-1
        exny(i,j)=vzny(i,j)*bym(i,j)-vym(i,j)*bzny(i,j)
      enddo
      enddo

      do j=2,jx-1
      do i=2,ix-1
        bzh(i,j)=bzh(i,j)
     &         -dt/dx(i)*(eynx(i,j)-eynx(i-1,j))
     &         +dt/dy(j)*(exny(i,j)-exny(i,j-1))
      enddo
      enddo

      return
      end
