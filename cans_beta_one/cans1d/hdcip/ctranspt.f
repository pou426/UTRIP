c======================================================================|
      subroutine ctranspt(byh,vxm,vyc,bxm,byc,dt,dx,ix)
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
c    vyc(ix): [double] velocity
c    byc(ix): [double] magnetic field
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

      dimension dx(ix)

      dimension byh(ix),ezm(ix)
      dimension vxm(ix),vyc(ix)
      dimension bxm(ix),byc(ix)
c----------------------------------------------------------------------|


      do i=1,ix-1
        ezm(i)=-(vxm(i)*byc(i)-vyc(i)*bxm(i))
      enddo

      do i=2,ix-1
        byh(i)=byh(i)+dt/dx(i)*(ezm(i)-ezm(i-1))
      enddo 


      return
      end
