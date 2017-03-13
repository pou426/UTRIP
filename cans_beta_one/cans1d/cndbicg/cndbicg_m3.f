c======================================================================|
      subroutine cndbicg_m3(ro,pr,mi,err,dt,gm,rkap0,margin
     &        ,bxm,by,bz,dx,dxm,ix)
c======================================================================|
c    
c NAME  cndbicg_m3
c 
c PURPOSE
c    solve heat conduction equation by SOR red & black method
c        * Spitzer type
c        * magnetic field 3 components
c 
c INPUTS & OUTPUTS
c    pr(ix): [double] pressure
c
c OUTPUTS
c    mi: [integer] number of iterations
c    err: [double] remained error for convergense
c
c INPUTS
c    ro(ix): [double] density
c    bxm(ix) : [double] magnetic field
c    by(ix) : [double] magnetic field
c    bz(ix) : [double] magnetic field
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    rkap0: [double] constant part of heat conduction coefficient
c    gm: [double] polytropic index gamma
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      parameter(eps=1.d-16)
      dimension ro(ix),pr(ix),te(ix)
      dimension by(ix),bz(ix)
      dimension bxm(ix)
      dimension cmat(ix,3),res(ix),src(ix)
      dimension dx(ix),dxm(ix)
      dimension work(ix,7),dmat(ix)
      integer r,rtld,p,phat,v,s,shat,t
c----------------------------------------------------------------------|

      call prtote(te,ro,pr,gm,ix)
      call ccfspt_m3(cmat,src,rkap0,gm,dt,te,ro,bxm,by,bz,dx,dxm,ix)
      call bndcmat(margin,cmat,ix)

      milim=ix
      eps0=sqrt(eps)

      mi=0
      call residue(anormf,src,margin,ix)
      call iludcmp(cmat,dmat,margin,ix)
      call bicgstab1(r,rtld,p,v,t,phat,shat,s,work,te,src,margin,ix)
 1000 continue
      mi=mi+1
      call bicgstab2(r,rtld,p,v,t,phat,shat,s,work,
     &     rho1,alpha,omega,te,res,cmat,dmat,margin,ix,mi)
      call residue(anorm,res,margin,ix)
      err=anorm/anormf
      if ((err.gt.eps0).and.(mi.lt.milim)) goto 1000

      call tetopr(ro,pr,te,gm,ix)

      return
      end
