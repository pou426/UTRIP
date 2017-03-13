c======================================================================|
      subroutine rtnewt_sh(ro,pr,vx,vy,gl
     &     ,de,ee,rx,ry,gm,mix,tolf,tolx,ix,i0,i1,jx,j0,j1)
c======================================================================|
      implicit double precision (a-h,o-z)

      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx),gl(ix,jx)
      dimension de(ix,jx),ee(ix,jx),rx(ix,jx),ry(ix,jx)
      dimension xx1(ix,jx),mflg(ix,jx)
      dimension rf2de(ix,jx),eede(ix,jx)

c----------------------------------------------------------------------|

      rk=gm/(gm-1.d0)

      do j=j0,j1
      do i=i0,i1

        mflg(i,j)=0

        rf2de(i,j)=(rx(i,j)**2+ry(i,j)**2)/(de(i,j)**2)
        eede(i,j)=ee(i,j)/de(i,j)

c--- Initial Guess

        glp=1.d0/sqrt(1.d0-vx(i,j)**2-vy(i,j)**2)
        xx1(i,j)=glp

      enddo
      enddo

c--- Iteration

      do mi=1,mix
      do j=j0,j1
      do i=i0,i1
      if (mflg(i,j).eq.0) then

        gg=xx1(i,j)

        psi=(eede(i,j)+1.d0-gg)/(rk*gg-1.d0/gg)
        hh=1.d0+rk*psi

        fv1=hh**2*(gg**2-1.d0)-rf2de(i,j)

        errf=abs(fv1)
        if (errf.le.tolf) mflg(i,j)=1

        c0=rk*gg-1.d0/gg
        dpsidgg= -1./c0*(1.d0+psi*(rk+1.d0/gg**2))

        a11=2.d0*hh*rk*(gg**2-1.d0)*dpsidgg+2.d0*hh**2*gg

        b11= 1./a11

        dxx1=-b11*fv1

        xx1(i,j)=xx1(i,j)+dxx1

        errx=abs(dxx1)
        if (errx.le.tolx) mflg(i,j)=1

      endif
      enddo
      enddo
      enddo

      do j=j0,j1
      do i=i0,i1
        gl(i,j)=xx1(i,j)
        ro(i,j)=de(i,j)/gl(i,j)
        psi=(eede(i,j)+1.d0-gl(i,j))/(rk*gl(i,j)-1.d0/gl(i,j))
        pr(i,j)=psi*ro(i,j)
        denom=(1.d0+rk*psi)*gl(i,j)*de(i,j)
        vx(i,j)=rx(i,j)/denom
        vy(i,j)=ry(i,j)/denom
      enddo
      enddo

      return
      end
