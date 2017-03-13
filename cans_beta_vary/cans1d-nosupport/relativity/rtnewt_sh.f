c======================================================================|
      subroutine rtnewt_sh(ro,pr,vx,gl
     &     ,de,ee,rx,gm,mix,tolf,tolx,ix,i0,i1)
c======================================================================|
      implicit double precision (a-h,o-z)

      dimension ro(ix),pr(ix)
      dimension vx(ix),gl(ix)
      dimension de(ix),ee(ix),rx(ix)
      dimension xx1(ix),mflg(ix)
      dimension rf2de(ix),eede(ix)

c----------------------------------------------------------------------|

      rk=gm/(gm-1.d0)

      do i=i0,i1

        mflg(i)=0

        rf2de(i)=(rx(i)**2)/(de(i)**2)
        eede(i)=ee(i)/de(i)

c--- Initial Guess

        glp=1.d0/sqrt(1.d0-vx(i)**2)
        xx1(i)=glp

      enddo

c--- Iteration

      do mi=1,mix
      do i=i0,i1
      if (mflg(i).eq.0) then

        gg=xx1(i)

        psi=(eede(i)+1.d0-gg)/(rk*gg-1.d0/gg)
        hh=1.d0+rk*psi

        fv1=hh**2*(gg**2-1.d0)-rf2de(i)

        errf=abs(fv1)
        if (errf.le.tolf) mflg(i)=1

        c0=rk*gg-1.d0/gg
        dpsidgg= -1./c0*(1.d0+psi*(rk+1.d0/gg**2))

        a11=2.d0*hh*rk*(gg**2-1.d0)*dpsidgg+2.d0*hh**2*gg

        b11= 1./a11

        dxx1=-b11*fv1

        xx1(i)=xx1(i)+dxx1

        errx=abs(dxx1)
        if (errx.le.tolx) mflg(i)=1

      endif
      enddo
      enddo

      do i=i0,i1
        gl(i)=xx1(i)
        ro(i)=de(i)/gl(i)
        psi=(eede(i)+1.d0-gl(i))/(rk*gl(i)-1.d0/gl(i))
        pr(i)=psi*ro(i)
        denom=(1.d0+rk*psi)*gl(i)*de(i)
        vx(i)=rx(i)/denom
      enddo

      return
      end
