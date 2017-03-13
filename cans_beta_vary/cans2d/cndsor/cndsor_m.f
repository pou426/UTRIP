c======================================================================|
      subroutine cndsor_m(ro,pr,mi,err,dt,gm,rkap0,bx,by
     &                   ,margin,dx,dxm,ix,dy,dym,jx)
c======================================================================|
c 
c NAME  cndsor_m
c 
c PURPOSE
c    solve heat conduction equation by SOR red & black method
c        * Spitzer type
c        * magnetic field ( 2 components )
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
c    ro(ix,jx): [double] density
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    rkap0: [double] constant part of heat conduction coefficient
c    gm: [double] polytropic index gamma
c    bx(ix,jx) : [double] magnetic field
c    by(ix,jx) : [double] magnetic field
c 
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      parameter(eps=1.d-16)
      dimension ro(ix,jx),pr(ix,jx),te(ix,jx)
      dimension bx(ix,jx),by(ix,jx)
      dimension cmat(ix,jx,5),res(ix,jx),src(ix,jx)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
        res(i,j)=0.d0
        src(i,j)=0.d0
      enddo
      enddo

      call prtote(te,ro,pr,gm,ix,jx)
      call ccfspt_m(cmat,src,rkap0,gm,dt,te,ro,bx,by
     &              ,dx,dxm,ix,dy,dym,jx)
      call bndcmat(margin,cmat,ix,jx)

      omsor=1.9
      milim=ix*jx
      eps0=sqrt(eps)

      mi=0
      call gtnorm(sumf,src,margin,ix,jx)

      mrb0=1
      mrb1=3-mrb0
1000  continue
      mi=mi+1
      call sorbr(te,res,mrb0,omsor,cmat,src,margin,ix,jx)
      call sorbr(te,res,mrb1,omsor,cmat,src,margin,ix,jx)
      call gtnorm(sum,res,margin,ix,jx)
      err=sqrt(sum/sumf)
      if ((err.gt.eps0).and.(mi.lt.milim)) goto 1000

      call tetopr(ro,pr,te,gm,ix,jx)

      return
      end
