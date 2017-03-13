c======================================================================|
      subroutine p_cndbicg(ro,pr,mi,err,dt,gm,rkap0
     &                    ,margin,dx,dxm,ix,dy,dym,jx
     &         ,igx,jgx,ipe,jpe,ipex,jpex,mper)
c======================================================================|
c    
c NAME  cndbicg
c    
c PURPOSE
c    solve heat conduction equation by biCG method
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
c    ro(ix,jx): [double] density
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    rkap0: [double] constant part of heat conduction coefficient
c    gm: [double] polytropic index gamma
c    ix,jx: [integer] dimension size
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      parameter(eps=1.d-16)
      dimension ro(ix,jx),pr(ix,jx),te(ix,jx)
      dimension cmat(ix,jx,5),res(ix,jx),src(ix,jx)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension work(ix,jx,7),dmat(ix,jx)
      integer   r,rtld,p,phat,v,s,shat,t

      include "mpif.h"
c----------------------------------------------------------------------|

      call prtote(te,ro,pr,gm,ix,jx)
      call ccfspt(cmat,src,rkap0,gm,dt,te,ro,dx,dxm,ix,dy,dym,jx)
c     call ccfunf(cmat,src,rkap0,gm,dt,te,ro,dx,dxm,ix,dy,dym,jx)
      call bndcmat(margin,cmat,ix,jx
     &           ,ipe,jpe,ipex,jpex)

      milim=igx*jgx
      eps0=sqrt(eps)

      mi=0
      call gtnorm(sumf,src,margin,ix,jx)
      call mpi_allreduce(sumf,sumfg,1,mpi_double_precision,mpi_sum
     &                      ,mpi_comm_world,merr)
      sumf=sumfg

      call iludcmp(cmat,dmat,margin,ix,jx)
      call bicgstab1(r,rtld,p,v,t,phat,shat,s,work,te,src
     &  ,margin,ix,jx)
1000  continue
      mi=mi+1
      call p_bicgstab2(r,rtld,p,v,t,phat,shat,s,work,
     &     rho1,alpha,omega,te,res,cmat,dmat,margin,ix,jx,mi,
     &     ipe,jpe,ipex,jpex,mper)
      call gtnorm(sum,res,margin,ix,jx)
      call mpi_allreduce(sum,sumg,1,mpi_double_precision,mpi_sum
     &                      ,mpi_comm_world,merr)
      sum=sumg
      err=sqrt(sum/sumf)

      if ((err.gt.eps0).and.(mi.lt.milim)) goto 1000

      call tetopr(ro,pr,te,gm,ix,jx)

      return
      end
