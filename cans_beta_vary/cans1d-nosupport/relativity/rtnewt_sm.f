c======================================================================|
      subroutine rtnewt_sm(ro,pr,vx,vy,gl
     &     ,de,ee,rx,ry,by,bx,gm,mix,tolf,tolx,ix,i0,i1)
c======================================================================|
      implicit double precision (a-h,o-z)

      dimension ro(ix),pr(ix)
      dimension vx(ix),vy(ix),gl(ix)
      dimension de(ix),ee(ix),rx(ix),ry(ix),by(ix)
      dimension bx(ix)
      dimension xx1(ix),xx2(ix),mflg(ix)
      dimension psim(ix),rhde(ix),rsde(ix),eede(ix)

c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi8i=1./pi/8.
      pi4i=1./pi/4.
c----------------------------------------------------------------------|

      rk=gm/(gm-1.d0)

      do i=i0,i1

        mflg(i)=0

        bf2=bx(i)**2+by(i)**2
        bf=sqrt(bf2)
        psim(i)=bf2*pi8i/de(i)

        rf2=rx(i)**2+ry(i)**2
        rh=(bx(i)/bf)*rx(i)+(by(i)/bf)*ry(i)
        rs2=rf2-rh**2
        if (rs2.lt.0.d0) rs2=0.d0
        rs=sqrt(rs2)

        rhde(i)=rh/de(i)
        rsde(i)=rs/de(i)
        eede(i)=ee(i)/de(i)

c--- Initial Guess

        vxp=vx(i)
        vyp=vy(i)
        vhp=(bx(i)*vxp+by(i)*vyp)/bf
        vfp2=vxp**2+vyp**2
        vsp2=vfp2-vhp**2
        if (vsp2.lt.0.d0) vsp2=0.d0
        vsp=sqrt(vsp2)

        vfp2=vsp2+vhp**2
        if (vfp2.gt.1.d0) vfp2=0.9999999999999999999999d0
        glp=1.d0/sqrt(1.d0-vfp2)

        xx1(i)=vsp
        xx2(i)=glp

      enddo

c--- Iteration

      do mi=1,mix
      do i=i0,i1
      if (mflg(i).eq.0) then

        vs=xx1(i)
        gg=xx2(i)

        psi=(eede(i)+1.d0-gg-psim(i)*(vs**2+1.d0))/(rk*gg-1.d0/gg)
        hh=1.d0+rk*psi

        fv1=hh*gg*vs+2.d0*psim(i)*vs-rsde(i)
        fv2=hh**2*(gg**2-1.d0)-(rsde(i)-2.d0*psim(i)*vs)**2-rhde(i)**2

        errf=abs(fv1)+abs(fv2)
        if (errf.le.tolf) mflg(i)=1

        c0=rk*gg-1.d0/gg
        dpsidvs = -2.d0*psim(i)*vs/c0
        dpsidgg= -1./c0*(1.d0+psi*(rk+1.d0/gg**2))

        a11=hh*gg+rk*gg*vs*dpsidvs +2.d0*psim(i)
        a12=hh*vs +rk*gg*vs*dpsidgg
        a21=2.d0*hh*rk*(gg**2-1.d0)*dpsidvs
     &        +4.d0*psim(i)*(rsde(i)-2.d0*psim(i)*vs)
        a22=2.d0*hh*rk*(gg**2-1.d0)*dpsidgg+2.d0*hh**2*gg

        det0= a11*a22-a12*a21

        b11= a22
        b21=-a12
        b12=-a21
        b22= a11

        dxx1=-(b11*fv1+b21*fv2)/det0
        dxx2=-(b12*fv1+b22*fv2)/det0

        xx1(i)=xx1(i)+dxx1
        xx2(i)=xx2(i)+dxx2

        errx=abs(dxx1)+abs(dxx2)
        if (errx.le.tolx) mflg(i)=1

      endif
      enddo
      enddo

      do i=i0,i1
        vs=xx1(i)
        gl(i)=xx2(i)
        ro(i)=de(i)/gl(i)
        psi=(eede(i)+1.d0-gl(i)-psim(i)*(vs**2+1.d0))
     &       /(rk*gl(i)-1.d0/gl(i))
        pr(i)=psi*ro(i)
        en = ro(i)+pr(i)/(gm-1.d0)
        vh=rhde(i)/(1.d0+rk*psi)/gl(i)
        bf2=bx(i)**2+by(i)**2
        bf=sqrt(bf2)
        denom=(en+pr(i))*gl(i)**2+bf2*pi4i
        vx(i) = (rx(i) + vh*bf*bx(i)*pi4i)/denom
        vy(i) = (ry(i) + vh*bf*by(i)*pi4i)/denom
      enddo

      return
      end
