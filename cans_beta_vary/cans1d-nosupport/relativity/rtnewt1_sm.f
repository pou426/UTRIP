c======================================================================|
      subroutine rtnewt1_sm(ro,pr,vx,vy,gl
     &     ,de,ee,rx,ry,by,bx,gm,mix,tolf,tolx,ix,i0,i1)
c======================================================================|
      implicit double precision (a-h,o-z)

      dimension ro(ix),pr(ix)
      dimension vx(ix),vy(ix),gl(ix)
      dimension de(ix),ee(ix),rx(ix),ry(ix),by(ix)
      dimension bx(ix)
      dimension xx1(ix),xx2(ix),xx3(ix),xx4(ix),mflg(ix)
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
        psip=pr(i)/ro(i)

        vfp2=vsp2+vhp**2
        if (vfp2.gt.1.d0) vfp2=0.9999999999999999999999d0
        glp=1.d0/sqrt(1.d0-vfp2)

        xx1(i)=vhp
        xx2(i)=vsp
        xx3(i)=glp
        xx4(i)=psip

      enddo

c--- Iteration

      do mi=1,mix
      do i=i0,i1
      if (mflg(i).eq.0) then

        vh=xx1(i)
        vs=xx2(i)
        gg=xx3(i)
        psi=xx4(i)

        hh=1.d0+rk*psi
        tt=1.d0-vh**2-vs**2

        fv1=hh*gg*vh-rhde(i)
        fv2=hh*gg*vs+2.d0*psim(i)*vs-rsde(i)
        fv3=hh*gg-psi/gg-1.d0 +psim(i)*(vs**2+1.d0)-eede(i)
        fv4=gg**2*tt-1.d0

        errf=abs(fv1)+abs(fv2)+abs(fv3)+abs(fv4)
        if (errf.le.tolf) mflg(i)=1

        c0=hh*gg+2.d0*psim(i)
        c1=1-rk*gg**2
        c2=(hh+rk*psi)*gg
        c3=hh*gg**2+psi

        det0= (hh**2*gg**2*tt*c1 + hh*gg*c2*(vh**2+vs**2)
     &     +hh*gg*tt*2.d0*psim(i)*(c1+rk*gg**2*vs**2) 
     &     +vh**2*2.d0*psim(i)*c2 )*2.

        b11=2*tt*c0*c1+2*vs**2*(c2+tt*2.d0*psim(i)*rk*gg**2)
        b21=-2*vs*vh*(c2+rk*gg**2*tt*2.d0*psim(i))
        b31=2*rk*gg**2*tt*vh*c0
        b41=-vh/gg**2*c0*c2

        b12=-2*vh*vs*c2
        b22=2*hh*gg*tt*c1+2.*vh**2*c2
        b32=2*rk*hh*gg**3*tt*vs
        b42=-hh*vs/gg*c2

        b13=2*hh*gg**2*vh*c1+2*gg*vh*2.d0*psim(i)*(c1+rk*gg**2*vs**2)
        b23=2*hh*gg**2*vs*c1-2*rk*gg**3*vh**2*vs*2.d0*psim(i)
        b33=2*rk*hh*gg**4*(vh**2+vs**2)+2*rk*gg**3*vh**2*2.d0*psim(i)
        b43=hh**2*gg*c1+hh*2.d0*psim(i)*(c1+rk*gg**2*vs**2)

        b14=2*hh*gg*vh*c3+2*vh*2.d0*psim(i)*(c3-hh*gg**2*vs**2)
        b24=2*hh*gg*vs*c3+2*hh*gg**2*vs*2.d0*psim(i)*(tt+vh**2)
        b34=-2*hh**2*gg**3*(tt+vh**2+vs**2)
     &         -2*hh*gg**2*2.d0*psim(i)*(tt+vh**2)
        b44=hh**2*c3+hh*2.d0*psim(i)/gg*(c3-hh*gg**2*vs**2)

        dxx1=-(b11*fv1+b21*fv2+b31*fv3+b41*fv4)/det0
        dxx2=-(b12*fv1+b22*fv2+b32*fv3+b42*fv4)/det0
        dxx3=-(b13*fv1+b23*fv2+b33*fv3+b43*fv4)/det0
        dxx4=-(b14*fv1+b24*fv2+b34*fv3+b44*fv4)/det0

        xx1(i)=xx1(i)+dxx1
        xx2(i)=xx2(i)+dxx2
        xx3(i)=xx3(i)+dxx3
        xx4(i)=xx4(i)+dxx4

        errx=abs(dxx1)+abs(dxx2)+abs(dxx3)+abs(dxx4)
        if (errx.le.tolx) mflg(i)=1

      endif
      enddo
      enddo


      do i=i0,i1
        vh=xx1(i)
        gl(i)=xx3(i)
        ro(i)=de(i)/gl(i)
        pr(i)=de(i)*xx4(i)/gl(i)
        en = ro(i)+pr(i)/(gm-1)
        bf2=bx(i)**2+by(i)**2
        bf=sqrt(bf2)
        denom=(en+pr(i))*gl(i)**2+bf2*pi4i
        vx(i) = (rx(i) + vh*bf*bx(i)*pi4i)/denom
        vy(i) = (ry(i) + vh*bf*by(i)*pi4i)/denom
      enddo

      return
      end
