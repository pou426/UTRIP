c======================================================================|
      subroutine rtnewt_sm(ro,pr,vx,vy,gl
     &     ,de,ee,rx,ry,bx,by,gm,mix,tolf,tolx
     &     ,ix,i0,i1,jx,j0,j1)
c======================================================================|
      implicit double precision (a-h,o-z)

      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx),gl(ix,jx)
      dimension de(ix,jx),ee(ix,jx),rx(ix,jx),ry(ix,jx),by(ix,jx)
      dimension bx(ix,jx)
      dimension xx1(ix,jx),xx2(ix,jx),mflg(ix,jx)
      dimension psim(ix,jx),rhde(ix,jx),rsde(ix,jx),eede(ix,jx)

c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi8i=1./pi/8.
      pi4i=1./pi/4.
c----------------------------------------------------------------------|

      rk=gm/(gm-1.d0)

      do j=j0,j1
      do i=i0,i1

        mflg(i,j)=0

        bf2=bx(i,j)**2+by(i,j)**2
        bf=sqrt(bf2)
        psim(i,j)=bf2*pi8i/de(i,j)

        rf2=rx(i,j)**2+ry(i,j)**2
        rh=(bx(i,j)/bf)*rx(i,j)+(by(i,j)/bf)*ry(i,j)
        rs2=rf2-rh**2
        if (rs2.lt.0.d0) rs2=0.d0
        rs=sqrt(rs2)

        rhde(i,j)=rh/de(i,j)
        rsde(i,j)=rs/de(i,j)
        eede(i,j)=ee(i,j)/de(i,j)

c--- Initial Guess

        vxp=vx(i,j)
        vyp=vy(i,j)
        vhp=(bx(i,j)*vxp+by(i,j)*vyp)/bf
        vfp2=vxp**2+vyp**2
        vsp2=vfp2-vhp**2
        if (vsp2.lt.0.d0) vsp2=0.d0
        vsp=sqrt(vsp2)

        vfp2=vsp2+vhp**2
        if (vfp2.gt.1.d0) vfp2=0.9999999999999999999999d0
        glp=1.d0/sqrt(1.d0-vfp2)

        xx1(i,j)=vsp
        xx2(i,j)=glp

      enddo
      enddo

c--- Iteration

      do mi=1,mix
      do j=i0,j1
      do i=i0,i1
      if (mflg(i,j).eq.0) then

        vs=xx1(i,j)
        gg=xx2(i,j)

        psi=(eede(i,j)+1.d0-gg-psim(i,j)*(vs**2+1.d0))/(rk*gg-1.d0/gg)
        hh=1.d0+rk*psi

        fv1=hh*gg*vs+2.d0*psim(i,j)*vs-rsde(i,j)
        fv2=hh**2*(gg**2-1.d0)-(rsde(i,j)-2.d0*psim(i,j)*vs)**2
     &       -rhde(i,j)**2

        errf=abs(fv1)+abs(fv2)
        if (errf.le.tolf) mflg(i,j)=1

        c0=rk*gg-1.d0/gg
        dpsidvs = -2.d0*psim(i,j)*vs/c0
        dpsidgg= -1./c0*(1.d0+psi*(rk+1.d0/gg**2))

        a11=hh*gg+rk*gg*vs*dpsidvs +2.d0*psim(i,j)
        a12=hh*vs +rk*gg*vs*dpsidgg
        a21=2.d0*hh*rk*(gg**2-1.d0)*dpsidvs
     &        +4.d0*psim(i,j)*(rsde(i,j)-2.d0*psim(i,j)*vs)
        a22=2.d0*hh*rk*(gg**2-1.d0)*dpsidgg+2.d0*hh**2*gg

        det0= a11*a22-a12*a21

        b11= a22
        b21=-a12
        b12=-a21
        b22= a11

        dxx1=-(b11*fv1+b21*fv2)/det0
        dxx2=-(b12*fv1+b22*fv2)/det0

        xx1(i,j)=xx1(i,j)+dxx1
        xx2(i,j)=xx2(i,j)+dxx2

        errx=abs(dxx1)+abs(dxx2)
        if (errx.le.tolx) mflg(i,j)=1

      endif
      enddo
      enddo
      enddo

      do j=j0,j1
      do i=i0,i1
        vs=xx1(i,j)
        gl(i,j)=xx2(i,j)
        ro(i,j)=de(i,j)/gl(i,j)
        psi=(eede(i,j)+1.d0-gl(i,j)-psim(i,j)*(vs**2+1.d0))
     &       /(rk*gl(i,j)-1.d0/gl(i,j))
        pr(i,j)=psi*ro(i,j)
        en = ro(i,j)+pr(i,j)/(gm-1.d0)
        vh=rhde(i,j)/(1.d0+rk*psi)/gl(i,j)
        bf2=bx(i,j)**2+by(i,j)**2
        bf=sqrt(bf2)
        denom=(en+pr(i,j))*gl(i,j)**2+bf2*pi4i
        vx(i,j) = (rx(i,j) + vh*bf*bx(i,j)*pi4i)/denom
        vy(i,j) = (ry(i,j) + vh*bf*by(i,j)*pi4i)/denom
      enddo
      enddo

      return
      end
