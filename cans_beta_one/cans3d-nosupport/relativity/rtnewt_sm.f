c======================================================================|
      subroutine rtnewt_sm(ro,pr,vx,vy,vz,gl
     &     ,de,ee,rx,ry,bx,by,bz,gm,mix,tolf,tolx
     &     ,ix,i0,i1,jx,j0,j1,kx,k0,k1)
c======================================================================|
      implicit double precision (a-h,o-z)

      dimension ro(ix,jx,kx),pr(ix,jx,kx),de(ix,jx,kx),ee(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx),gl(ix,jx,kx)
      dimension rx(ix,jx,kx),ry(ix,jx,kx),rz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension xx1(ix,jx,kx),xx2(ix,jx,kx),mflg(ix,jx,kx)
      dimension eede(ix,jx,kx),psim(ix,jx,kx)
      dimension rhde(ix,jx,kx),rsde(ix,jx,kx)

c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi8i=1./pi/8.
      pi4i=1./pi/4.
c----------------------------------------------------------------------|

      rk=gm/(gm-1.d0)

      do k=k0,k1
      do j=j0,j1
      do i=i0,i1

        mflg(i,j,k)=0

        bf2=bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2
        bf=sqrt(bf2)
        psim(i,j,k)=bf2*pi8i/de(i,j,k)

        rf2=rx(i,j,k)**2+ry(i,j,k)**2+rz(i,j,k)**2
        rh=(bx(i,j,k)/bf)*rx(i,j,k)
     &    +(by(i,j,k)/bf)*ry(i,j,k)
     &    +(bz(i,j,k)/bf)*rz(i,j,k)
        rs2=rf2-rh**2
        if (rs2.lt.0.d0) rs2=0.d0
        rs=sqrt(rs2)

        rhde(i,j,k)=rh/de(i,j,k)
        rsde(i,j,k)=rs/de(i,j,k)
        eede(i,j,k)=ee(i,j,k)/de(i,j,k)

c--- Initial Guess

        vxp=vx(i,j,k)
        vyp=vy(i,j,k)
        vzp=vz(i,j,k)
        vhp=(bx(i,j,k)*vxp+by(i,j,k)*vyp+bz(i,j,k)*vzp)/bf
        vfp2=vxp**2+vyp**2+vzp**2
        vsp2=vfp2-vhp**2
        if (vsp2.lt.0.d0) vsp2=0.d0
        vsp=sqrt(vsp2)

        vfp2=vsp2+vhp**2
        if (vfp2.gt.1.d0) vfp2=0.9999999999999999999999d0
        glp=1.d0/sqrt(1.d0-vfp2)

        xx1(i,j,k)=vsp
        xx2(i,j,k)=glp

      enddo
      enddo
      enddo

c--- Iteration

      do mi=1,mix
      do k=k0,k1
      do j=j0,j1
      do i=i0,i1
      if (mflg(i,j,k).eq.0) then

        vs=xx1(i,j,k)
        gg=xx2(i,j,k)

        psi=(eede(i,j,k)+1.d0-gg-psim(i,j,k)*(vs**2+1.d0))
     &       /(rk*gg-1.d0/gg)
        hh=1.d0+rk*psi

        fv1=hh*gg*vs+2.d0*psim(i,j,k)*vs-rsde(i,j,k)
        fv2=hh**2*(gg**2-1.d0)-(rsde(i,j,k)-2.d0*psim(i,j,k)*vs)**2
     &       -rhde(i,j,k)**2

        errf=abs(fv1)+abs(fv2)
        if (errf.le.tolf) mflg(i,j,k)=1

        c0=rk*gg-1.d0/gg
        dpsidvs = -2.d0*psim(i,j,k)*vs/c0
        dpsidgg= -1./c0*(1.d0+psi*(rk+1.d0/gg**2))

        a11=hh*gg+rk*gg*vs*dpsidvs +2.d0*psim(i,j,k)
        a12=hh*vs +rk*gg*vs*dpsidgg
        a21=2.d0*hh*rk*(gg**2-1.d0)*dpsidvs
     &        +4.d0*psim(i,j,k)*(rsde(i,j,k)-2.d0*psim(i,j,k)*vs)
        a22=2.d0*hh*rk*(gg**2-1.d0)*dpsidgg+2.d0*hh**2*gg

        det0= a11*a22-a12*a21

        b11= a22
        b21=-a12
        b12=-a21
        b22= a11

        dxx1=-(b11*fv1+b21*fv2)/det0
        dxx2=-(b12*fv1+b22*fv2)/det0

        xx1(i,j,k)=xx1(i,j,k)+dxx1
        xx2(i,j,k)=xx2(i,j,k)+dxx2

        errx=abs(dxx1)+abs(dxx2)
        if (errx.le.tolx) mflg(i,j,k)=1

      endif
      enddo
      enddo
      enddo
      enddo

      do k=k0,k1
      do j=j0,j1
      do i=i0,i1
        vs=xx1(i,j,k)
        gl(i,j,k)=xx2(i,j,k)
        ro(i,j,k)=de(i,j,k)/gl(i,j,k)
        psi=(eede(i,j,k)+1.d0-gl(i,j,k)-psim(i,j,k)*(vs**2+1.d0))
     &       /(rk*gl(i,j,k)-1.d0/gl(i,j,k))
        pr(i,j,k)=psi*ro(i,j,k)
        en = ro(i,j,k)+pr(i,j,k)/(gm-1.d0)
        vh=rhde(i,j,k)/(1.d0+rk*psi)/gl(i,j,k)
        bf2=bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2
        bf=sqrt(bf2)
        denom=(en+pr(i,j,k))*gl(i,j,k)**2+bf2*pi4i
        vx(i,j,k) = (rx(i,j,k) + vh*bf*bx(i,j,k)*pi4i)/denom
        vy(i,j,k) = (ry(i,j,k) + vh*bf*by(i,j,k)*pi4i)/denom
        vz(i,j,k) = (rz(i,j,k) + vh*bf*bz(i,j,k)*pi4i)/denom
      enddo
      enddo
      enddo

      return
      end
