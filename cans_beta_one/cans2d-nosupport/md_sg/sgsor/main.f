c======================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
      parameter (ng=6,nx=2**ng)
      parameter (ix=nx+3,jx=nx+3)
c     parameter (ix=10,jx=10)

      dimension x(ix),y(jx)
      dimension dx(ix),dy(jx)
      dimension dxm(ix),dym(jx)
      dimension ro(ix,jx),gx(ix,jx),gy(ix,jx)
      dimension gxm(ix,jx),gym(ix,jx)

      dimension uod(nx+2,nx+2)
      dimension gp(ix,jx)

c======================================================================|
c     prologue
c======================================================================|
      margin=2

      pi=4.0d0*datan(1.0d0)
      do i=1,ix
         x(i)=2.0d0*pi*float(i)/float(ix)
         dx(i)=2.d0*pi/float(ix)
         dxm(i)=2.d0*pi/float(ix)
      enddo
      do j=1,jx
         y(j)=2.0d0*pi*float(j)/float(jx)
         dy(j)=2.d0*pi/float(jx)
         dym(j)=2.d0*pi/float(jx)
      enddo
      do j=1,jx
      do i=1,ix
c        ro(i,j)=1.0d0+0.1d0*sin(x(i))*sin(y(j))
         w=0.5
         rr2=(x(i)-1.5*pi)**2+(y(j)-pi)**2
         ro(i,j)=0.1d0*exp(-rr2/w**2)
c        ro(i,j)=-sin(y(j))
      enddo
      enddo

c======================================================================|
c     main part
c======================================================================|
      call selfgsor(ro,gp,mi,err,margin,dx,dxm,ix,dy,dym,jx)
      call gptogg(gx,gy,gxm,gym,gp,dx,dxm,dy,dym,ix,jx)

c======================================================================|
c     epilogue
c======================================================================|

      mf_params=9
      call dacdefparam(mf_params,'params.txt')

      call dacputparamc(mf_params,'comment','cans2d md_sg')
      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'jx',jx)
      mf_x=11
      call dacdef1d(mf_x,'x.dac',6,ix)
      write(mf_x) x
      mf_y=12
      call dacdef1d(mf_y,'y.dac',6,jx)
      write(mf_y) y
      mf_ro=13
      call dacdef2d(mf_ro,'ro.dac',6,ix,jx)
      write(mf_ro) ro
      mf_gx=14
      call dacdef2d(mf_gx,'gx.dac',6,ix,jx)
      write(mf_gx) gx
      mf_gy=15
      call dacdef2d(mf_gy,'gy.dac',6,ix,jx)
      write(mf_gy) gy
      mf_gp=16
      call dacdef2d(mf_gp,'gp.dac',6,ix,jx)
      write(mf_gp) gp
      mf_uod=17
      call dacdef2d(mf_uod,'uod.dac',6,nx+2,nx+2)
      write(mf_uod) uod


      stop
      end
