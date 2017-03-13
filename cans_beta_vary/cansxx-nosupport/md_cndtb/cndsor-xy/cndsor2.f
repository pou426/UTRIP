c======================================================================|
      subroutine cndsor2(te,err,eps0,src,cmat,omsor
     &      ,mi,milim,ix,jx,kx,mfdim,margar)
c======================================================================|
c    
c NAME  cndbicg
c    
c PURPOSE
c    solve heat conduction equation by biCG method
c        * Spitzer type
c     
c INPUTS & OUTPUTS
c    pr(ix,jx,kx): [double] pressure
c
c OUTPUTS
c    mi: [integer] number of iterations
c    err: [double] remained error for convergense
c 
c INPUTS
c    ro(ix,jx,kx): [double] density
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    rkap0: [double] constant part of heat conduction coefficient
c    gm: [double] polytropic index gamma
c    ix,jx,kx: [integer] dimension size
c    
c HISTORY
c    written 2004-3-26 K. Nakamura
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension mfdim(3),margar(3)
      dimension te(ix,jx,kx)
      dimension cmat(ix,jx,kx,7),res(ix,jx,kx),src(ix,jx,kx)
      dimension work(ix,jx,kx,7)
c----------------------------------------------------------------------|

      mi=0
      call gtnorm2(sumf,src,margar,ix,jx,kx)

      mrb0=0
      mrb1=1
1000  continue
      mi=mi+1
      call sorbr2(te,res,mrb0,omsor,cmat,src,ix,jx,kx,mfdim,margar)
      call sorbr2(te,res,mrb1,omsor,cmat,src,ix,jx,kx,mfdim,margar)
      call gtnorm2(sum,res,margar,ix,jx,kx)
      err=sqrt(sum/sumf)
      if ((err.gt.eps0).and.(mi.lt.milim)) goto 1000


      return
      end
