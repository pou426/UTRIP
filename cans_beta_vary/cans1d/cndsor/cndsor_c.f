c======================================================================|
      subroutine cndsor_c(ro,pr,mi,err,dt,gm,rkap0,margin
     &    ,sc,scm,dx,dxm,ix)
c======================================================================|
c
c NAME  cndsor_c
c
c PURPOSE
c    solve heat conduction equation by SOR red & black method
c        * Spitzer type
c        * non-uniform cross section
c
c INPUTS & OUTPUTS
c    pr(ix): [double] pressure
c    ro(ix): [double] density
c
c OUTPUTS
c    mi: [integer] number of iterations
c    err: [double] remained error for convergense
c
c INPUTS
c    ix: [integer] dimension size
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    rkap0: [double] constant part of heat conduction coefficient
c    gm: [double] polytropic index gamma
c    sc(ix), scm(ix) : [double] cross section
c    margin: [integer] margin, i.e. # of grid points outside the boundary
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      parameter(eps=1.d-16)
      dimension ro(ix),pr(ix),te(ix)
      dimension cmat(ix,3),res(ix),src(ix)
      dimension dx(ix),dxm(ix)
      dimension sc(ix),scm(ix)
c----------------------------------------------------------------------|

      call prtote(te,ro,pr,gm,ix)
      call ccfspt_c(cmat,src,rkap0,gm,dt,te,ro,sc,scm,dx,dxm,ix)
c     call ccfunf_c(cmat,src,rkap0,gm,dt,te,ro,sc,scm,dx,dxm,ix)
      call bndcmat(margin,cmat,ix)

      omsor=1.9
      milim=ix*ix
      eps0=sqrt(eps)

      mi=0
      call residue(anormf,src,margin,ix)
1000  continue
      mi=mi+1
      call sorbr(te,res,1,omsor,cmat,src,margin,ix)
      call sorbr(te,res,2,omsor,cmat,src,margin,ix)
      call residue(anorm,res,margin,ix)
      err=anorm/anormf
      if ((err.gt.eps0).and.(mi.lt.milim)) goto 1000

      call tetopr(ro,pr,te,gm,ix)

      return
      end
