c======================================================================|
      subroutine p_cndbicg(ro,pr,mi,err,dt,gm,rkap0
     &     ,margin,dx,dxm,ix,dy,dym,jx,dz,dzm,kx
     &         ,igx,jgx,kgx,ipe,jpe,kpe,ipex,jpex,kpex,mper)
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
      parameter(eps=1.d-16)
      dimension ro(ix,jx,kx),pr(ix,jx,kx),te(ix,jx,kx)
      dimension cmat(ix,jx,kx,7),res(ix,jx,kx),src(ix,jx,kx)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension dz(kx),dzm(kx)
      dimension work(ix,jx,kx,7),dmat(ix,jx,kx)
      integer   r,rtld,p,phat,v,s,shat,t
      dimension mper(3)
      
      include "mpif.h"
c----------------------------------------------------------------------|

      call prtote(te,ro,pr,gm,ix,jx,kx)
      call ccfspt(cmat,src,rkap0,gm,dt,te,ro,dx,dxm,ix,dy,dym,jx
     &     ,dz,dzm,kx)
c     call ccfunf(cmat,src,rkap0,gm,dt,te,ro,dx,dxm,ix,dy,dym,jx
c     &     ,dz,dzm,kx)
      call bndcmat(margin,cmat,ix,jx,kx
     &           ,ipe,jpe,kpe,ipex,jpex,kpex)

      milim=igx*jgx*kgx
      eps0=sqrt(eps)

      mi=0
      call gtnorm(sumf,src,margin,ix,jx,kx)
      call mpi_allreduce(sumf,sumfg,1,mpi_double_precision,
     &    mpi_sum,mpi_comm_world,merr)
      sumf=sumfg

      call iludcmp(cmat,dmat,margin,ix,jx,kx)
      call bicgstab1(r,rtld,p,v,t,phat,shat,s,work,te,src,margin,ix,jx,
     &     kx)
1000  continue
      mi=mi+1
      call p_bicgstab2(r,rtld,p,v,t,phat,shat,s,work,
     &     rho1,alpha,omega,te,res,cmat,dmat,margin,ix,jx,kx,mi,
     &     ipe,jpe,kpe,ipex,jpex,kpex,mper)
      call gtnorm(sum,res,margin,ix,jx,kx)
      call mpi_allreduce(sum,sumg,1,mpi_double_precision,mpi_sum,
     &    mpi_comm_world,merr)
      sum=sumg
      err=sqrt(sum/sumf)

      if ((err.gt.eps0).and.(mi.lt.milim)) goto 1000

      call tetopr(ro,pr,te,gm,ix,jx,kx)

      return
      end
