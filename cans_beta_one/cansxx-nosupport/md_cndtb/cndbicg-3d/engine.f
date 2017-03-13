c======================================================================|
      subroutine engine(ro,pr,mi,err,dt,gm,rkap0
     &                    ,margin,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)
c======================================================================|
      implicit double precision (a-h,o-z)
      parameter(eps=1.d-16)
      dimension ro(ix,jx,kx),pr(ix,jx,kx),te(ix,jx,kx)
      dimension cmat(ix,jx,kx,7),res(ix,jx,kx),src(ix,jx,kx)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension dz(kx),dzm(kx)
      dimension work(ix,jx,kx,7),dmat(ix,jx,kx)
      integer   r,rtld,p,phat,v,s,shat,t
c----------------------------------------------------------------------|

      call eh2pr(pr,eh,gm,ix,jx,kx)
      call prtote(te,ro,pr,gm,ix,jx,kx)
      call ccfspt(cmat,src,rkap0,gm,dt,te,ro,dx,dxm,ix,dy,dym,jx
     &     ,dz,dzm,kx)
c     call ccfunf(cmat,src,rkap0,gm,dt,te,ro,dx,dxm,ix,dy,dym,jx
c     &     ,dz,dzm,kx)
      call bndcmat(margin,cmat,ix,jx,kx)

      milim=ix*jx*kx
      eps0=sqrt(eps)

      call cndbicg2(te,mi,err,eps0,src,cmat
     &                    ,margin,dx,dxm,ix,dy,dym,jx,dz,dzm,kx,milim)

c      mi=0
c      call gtnorm(sumf,src,margin,ix,jx,kx)

c      call iludcmp(cmat,dmat,margin,ix,jx,kx)
c      call bicgstab1(r,rtld,p,v,t,phat,shat,s,work,te,src,margin,ix,jx,
c     &     kx)
c1000  continue
c      mi=mi+1
c      call bicgstab2(r,rtld,p,v,t,phat,shat,s,work,
c     &     rho1,alpha,omega,te,res,cmat,dmat,margin,ix,jx,kx,mi)
c      call gtnorm(sum,res,margin,ix,jx,kx)
c      err=sqrt(sum/sumf)

c      if ((err.gt.eps0).and.(mi.lt.milim)) goto 1000

      call tetopr(ro,pr,te,gm,ix,jx,kx)

      return
      end
