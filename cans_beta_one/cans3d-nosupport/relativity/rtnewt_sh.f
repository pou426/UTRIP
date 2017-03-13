c======================================================================|
      subroutine rtnewt_sh(ro,pr,vx,vy,vz,gl
     &     ,de,ee,rx,ry,rz,gm,mix,tolf,tolx
     &     ,ix,i0,i1,jx,j0,j1,kx,k0,k1)
c======================================================================|
      implicit double precision (a-h,o-z)

      dimension ro(ix,jx,kx),pr(ix,jx,kx),de(ix,jx,kx),ee(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx),gl(ix,jx,kx)
      dimension rx(ix,jx,kx),ry(ix,jx,kx),rz(ix,jx,kx)
      dimension xx1(ix,jx,kx),mflg(ix,jx,kx)
      dimension rf2de(ix,jx,kx),eede(ix,jx,kx)

c----------------------------------------------------------------------|

      rk=gm/(gm-1.d0)

      do k=k0,k1
      do j=j0,j1
      do i=i0,i1

        mflg(i,j,k)=0

        rf2=rx(i,j,k)**2+ry(i,j,k)**2+rz(i,j,k)**2
        rf2de(i,j,k)=rf2/(de(i,j,k)**2)
        eede(i,j,k)=ee(i,j,k)/de(i,j,k)

c--- Initial Guess

        vf2=vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
        glp=1.d0/sqrt(1.d0-vf2)
        xx1(i,j,k)=glp

      enddo
      enddo
      enddo

c--- Iteration

      do mi=1,mix
      do k=k0,k1
      do j=j0,j1
      do i=i0,i1
      if (mflg(i,j,k).eq.0) then

        gg=xx1(i,j,k)

        psi=(eede(i,j,k)+1.d0-gg)/(rk*gg-1.d0/gg)
        hh=1.d0+rk*psi

        fv1=hh**2*(gg**2-1.d0)-rf2de(i,j,k)

        errf=abs(fv1)
        if (errf.le.tolf) mflg(i,j,k)=1

        c0=rk*gg-1.d0/gg
        dpsidgg= -1./c0*(1.d0+psi*(rk+1.d0/gg**2))

        a11=2.d0*hh*rk*(gg**2-1.d0)*dpsidgg+2.d0*hh**2*gg

        b11= 1./a11

        dxx1=-b11*fv1

        xx1(i,j,k)=xx1(i,j,k)+dxx1

        errx=abs(dxx1)
        if (errx.le.tolx) mflg(i,j,k)=1

      endif
      enddo
      enddo
      enddo
      enddo

      do k=k0,k1
      do j=j0,j1
      do i=i0,i1
        gl(i,j,k)=xx1(i,j,k)
        ro(i,j,k)=de(i,j,k)/gl(i,j,k)
        psi=(eede(i,j,k)+1.d0-gl(i,j,k))/(rk*gl(i,j,k)-1.d0/gl(i,j,k))
        pr(i,j,k)=psi*ro(i,j,k)
        denom=(1.d0+rk*psi)*gl(i,j,k)*de(i,j,k)
        vx(i,j,k)=rx(i,j,k)/denom
        vy(i,j,k)=ry(i,j,k)/denom
        vz(i,j,k)=rz(i,j,k)/denom
      enddo
      enddo
      enddo

      return
      end
