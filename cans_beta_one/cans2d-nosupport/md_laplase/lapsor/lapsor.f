c======================================================================|
      subroutine lapsor(pot,x,y,mi,err
     &          ,margin,dx,dxm,ix,dy,dym,jx)
c======================================================================|
c
c NAME  cndsor
c
c PURPOSE
c    solve heat conduction equation by SOR red & black method
c        * Spitzer type
c
c INPUTS & OUTPUTS
c    pr(ix,jx): [double] pressure
c                
c OUTPUTS
c    mi: [integer] number of iterations
c    err: [double] remained error for convergense
c
c INPUTS
c    ix,jx: [integer] dimension size
c    dx(ix), dxm(ix): [double] grid spacing
c     
c HISTORY
c    written 2002-3-1 T. Yokoyama
c     
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      parameter(eps=1.d-16)
      dimension pot(ix,jx)
      dimension cmat(ix,jx,5),res(ix,jx),src(ix,jx)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension x(ix),y(jx)
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
        res(i,j)=0.d0       
        src(i,j)=0.d0       
      enddo
      enddo

      call ccflap(cmat,src,dx,dxm,ix,dy,dym,jx)
      call bnd(pot,x,y,margin,ix,jx)

      omsor=1.9
      milim=ix*jx
      eps0=eps

      mi=0
      call gtnorm(sumf,src,margin,ix,jx)

      mrb0=1
      mrb1=3-mrb0
1000  continue
      mi=mi+1
      call sorbr(pot,res,mrb0,omsor,cmat,src,margin,ix,jx)
      call sorbr(pot,res,mrb1,omsor,cmat,src,margin,ix,jx)
      call bnd(pot,x,y,margin,ix,jx)

      call gtnorm(sum,res,margin,ix,jx)
c     err=sqrt(sum/sumf)
      err=sum
      if ((err.gt.eps0).and.(mi.lt.milim)) goto 1000

      return
      end
