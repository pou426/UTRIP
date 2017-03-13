c======================================================================|
      subroutine selfgbicg(ro,gp,mi,err,margin,dx,dxm,ix,dy,dym,jx)
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
      dimension ro(ix,jx),gp(ix,jx)
      dimension rosrc(ix,jx)
      dimension cmat(ix,jx,5),res(ix,jx),src(ix,jx)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension work(ix,jx,7),dmat(ix,jx)
      integer   r,rtld,p,phat,v,s,shat,t
c----------------------------------------------------------------------|

      rmass=0.
      vol=0.
      do i=1,ix
      do j=1,jx
        rmass=rmass+ro(i,j)*dx(i)*dy(j)
        vol=vol+dx(i)*dy(j)
      enddo
      enddo
      roavg=rmass/vol

      do i=1,ix
      do j=1,jx
        rosrc(i,j)=ro(i,j)-roavg
      enddo
      enddo

      call psncf(cmat,src,rosrc,dx,dxm,ix,dy,dym,jx)
      call bndcmat(margin,cmat,ix,jx)

      milim=ix*jx
      eps0=sqrt(eps)

      mi=0
      call gtnorm(sumf,src,margin,ix,jx)

      call iludcmp(cmat,dmat,margin,ix,jx)
      call bicgstab1(r,rtld,p,v,t,phat,shat,s,work,gp,src,margin,ix,jx)
1000  continue
      mi=mi+1
      call bicgstab2(r,rtld,p,v,t,phat,shat,s,work,
     &     rho1,alpha,omega,gp,res,cmat,dmat,margin,ix,jx,mi)

      call gtnorm(sum,res,margin,ix,jx)
      err=sqrt(sum/sumf)

      if ((err.gt.eps0).and.(mi.lt.milim)) goto 1000

      return
      end
