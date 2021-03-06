c======================================================================|
      subroutine selfgsor(ro,gp,mi,err,margin,dx,dxm,ix,dy,dym,jx)
c======================================================================|
c
c NAME  sgsor
c
c PURPOSE
c    solve heat conduction equation by SOR red & black method
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
c    ix,jx: [integer] dimension size
c    ro(ix,jx): [double] density
c    gm: [double] polytropic index gamma
c    rkap0: [double] constant part of heat conduction coefficient
c    dx(ix), dxm(ix): [double] grid spacing
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
      dimension dx(ix),dxm(ix),dy(jx),dym(jx)
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix
        res(i,j)=0.d0       
        src(i,j)=0.d0       
        gp(i,j)=0.d0       
      enddo
      enddo
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

      omsor=1.8
      milim=ix*jx
      eps0=sqrt(eps)

      mi=0
      call gtnorm(sumf,src,margin,ix,jx)

      mrb0=1
      mrb1=3-mrb0
1000  continue
      mi=mi+1
      call sorbr1(gp,res,mrb0,omsor,cmat,src,margin,ix,jx,mi)
      call sorbr1(gp,res,mrb1,omsor,cmat,src,margin,ix,jx,mi)

      call bdperx(margin,margin,gp,ix,jx)
      call bdpery(margin,margin,gp,ix,jx)

      call gtnorm(sum,res,margin,ix,jx)
      err=sqrt(sum/sumf)
      if ((err.gt.eps0).and.(mi.lt.milim)) goto 1000

      return
      end
