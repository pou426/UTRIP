c======================================================================|
      subroutine cndsor_m(ro,pr,mi,err,dt,gm,rkap0,bx,by,bz
     &         ,margin,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)
c======================================================================|
c 
c NAME  cndsor_m
c 
c PURPOSE
c    solve heat conduction equation by SOR red & black method
c        * Spitzer type
c        * magnetic field
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
c    dt: [double] delta time 
c    rkap0: [double] constant part of heat conduction coefficient
c    gm: [double] polytropic index gamma
c    bx(ix,jx,kx) : [double] magnetic field
c    by(ix,jx,kx) : [double] magnetic field 
c     
c HISTORY
c    written 2002-3-1 T. Yokoyama
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
      call ccfspt_m(cmat,src,rkap0,gm,dt,te,ro,bx,by,bz
     &            ,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)
      call bndcmat(margin,cmat,ix,jx,kx)

      omsor=1.9
      milim=ix*jx*kx
      eps0=sqrt(eps)

      mi=0
      call gtnorm(sumf,src,margin,ix,jx,kx)

      mrb0=1
      mrb1=3-mrb0
1000  continue
      mi=mi+1
      call sorbr(te,res,mrb0,omsor,cmat,src,margin,ix,jx,kx)
      call sorbr(te,res,mrb1,omsor,cmat,src,margin,ix,jx,kx)
      call gtnorm(sum,res,margin,ix,jx,kx)
      err=sqrt(sum/sumf)
      if ((err.gt.eps0).and.(mi.lt.milim)) goto 1000

      call tetopr(ro,pr,te,gm,ix,jx,kx)

      return
      end
