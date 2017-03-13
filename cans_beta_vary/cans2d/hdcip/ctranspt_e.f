c======================================================================|
      subroutine ctranspt_e(bxmh,bymh,vxc,vyc,bxc,byc,etm,czc
     &   ,dt,dx,dy,ix,jx)
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

      dimension bxmh(ix,jx),bymh(ix,jx),ezc(ix,jx)
      dimension vxc(ix,jx),vyc(ix,jx)
      dimension bxc(ix,jx),byc(ix,jx)
      dimension czc(ix,jx),etm(ix,jx)
c----------------------------------------------------------------------|


      do j=1,jx-1
      do i=1,ix-1
        ezc(i,j)=-(vxc(i,j)*byc(i,j)-vyc(i,j)*bxc(i,j))
     &           +etm(i,j)*czc(i,j)
      enddo
      enddo

      do j=1,jx-1
      do i=2,ix-1
        bymh(i,j)=bymh(i,j)+dt/dx(i)*(ezc(i,j)-ezc(i-1,j))
      enddo
      enddo

      do j=2,jx-1
      do i=1,ix-1
        bxmh(i,j)=bxmh(i,j)-dt/dy(j)*(ezc(i,j)-ezc(i,j-1))
      enddo
      enddo


      return
      end
