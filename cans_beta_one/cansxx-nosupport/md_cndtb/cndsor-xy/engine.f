c======================================================================|
      subroutine engine(ro,pr,err,dt,gm,rkap0
     &       ,dx,dxm,dy,dym,dz,dzm,mi,ix,jx,kx,mfdim,margar)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension mfdim(3),margar(3)
      parameter(eps=1.d-16)
      dimension ro(ix,jx,kx),pr(ix,jx,kx),te(ix,jx,kx)
      dimension cmat(ix,jx,kx,7),res(ix,jx,kx),src(ix,jx,kx)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension dz(kx),dzm(kx)
      dimension work(ix,jx,kx,7),dmat(ix,jx,kx)
      integer   r,rtld,p,phat,v,s,shat,t
c----------------------------------------------------------------------|

      call prtote(te,ro,pr,gm,ix,jx,kx)
      call ccfspt(cmat,src,rkap0,gm,dt,te,ro,dx,dxm,ix,dy,dym,jx
     &     ,dz,dzm,kx)
      call bndcmat(cmat,ix,jx,kx,mfdim,margar)

      omsor=1.9d0
      milim=ix*jx*kx
      eps0=sqrt(eps)

      call cndsor2(te,err,eps0,src,cmat,omsor
     &      ,mi,milim,ix,jx,kx,mfdim,margar)

      call tetopr(ro,pr,te,gm,ix,jx,kx)

      return
      end
