c======================================================================|
      subroutine cndbicg_mc(ro,pr,mi,err,dt,gm,rkap0,bx,by,bz
     &                    ,margin,x,xh,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)
c======================================================================|
c    
c NAME  cndbicg_mc
c    
c PURPOSE
c    solve heat conduction equation by biCG method
c        * Spitzer type
c        * magnetic field
c        * Cylindrical coordinate
c     
c INPUTS & OUTPUTS
c    pr(ix,jx,kx): [double] pressure
c
c OUTPUTS
c    mi: [integer] number of iterations
c    err: [double] remained error for convergense
c 
c INPUTS
c    ix,jx,kx: [integer] dimension size
c    ro(ix,jx,kx): [double] density
c    dx(ix), dxm(ix): [double] grid spacing
c    x(ix),xh(ix): [double] coordinate
c    dt: [double] delta time
c    rkap0: [double] constant part of heat conduction coefficient
c    gm: [double] polytropic index gamma
c    bx(ix,jx,kx) : [double] magnetic field
c    by(ix,jx,kx) : [double] magnetic field 
c    bz(ix,jx,kx) : [double] magnetic field 
c    margin: [integer] margin, i.e. # of grid points outside the boundary
c    
c HISTORY
c    written 2004-3-26 K. Nakamura
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      parameter(eps=1.d-16)
      dimension ro(ix,jx,kx),pr(ix,jx,kx),te(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension cmat(ix,jx,kx,7),res(ix,jx,kx),src(ix,jx,kx)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension dz(kx),dzm(kx)
      dimension x(ix),xh(ix)
      dimension work(ix,jx,kx,7),dmat(ix,jx,kx)
      integer   r,rtld,p,phat,v,s,shat,t
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
        res(i,j,k)=0.d0
        src(i,j,k)=0.d0
      enddo
      enddo
      enddo

      call prtote(te,ro,pr,gm,ix,jx,kx)
      call ccfspt_mc(cmat,src,rkap0,gm,dt,te,ro,bx,by,bz
     &              ,x,xh,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)
      call bndcmat(margin,cmat,ix,jx,kx)

      milim=ix*jx*kx
      eps0=sqrt(eps)

      mi=0
      call gtnorm(sumf,src,margin,ix,jx,kx)

      call iludcmp(cmat,dmat,margin,ix,jx,kx)
      call bicgstab1(r,rtld,p,v,t,phat,shat,s,work,
     &               te,src,margin,ix,jx,kx)
1000  continue
      mi=mi+1
      call bicgstab2(r,rtld,p,v,t,phat,shat,s,work,
     &               rho1,alpha,omega,te,res,cmat,dmat,margin,
     &               ix,jx,kx,mi)
      call gtnorm(sum,res,margin,ix,jx,kx)
      err=sqrt(sum/sumf)
      if ((err.gt.eps0).and.(mi.lt.milim)) goto 1000

      call tetopr(ro,pr,te,gm,ix,jx,kx)

      return
      end
