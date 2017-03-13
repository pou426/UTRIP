c======================================================================|
      subroutine rtbis_sh(ro,pr,vx,gl,de,ee,rx,gm,mix,tolf,tolx,xxmax
     &     ,ix,i0,i1)
c======================================================================|
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension de(ix),ee(ix),rx(ix)
      dimension ro(ix),pr(ix),vx(ix),gl(ix)
      dimension eede(ix),rf2de(ix)
      dimension fmid(ix),xx(ix),dxx(ix),mflg(ix)
c----------------------------------------------------------------------|
      func1(gg,eede0,rk)=(eede0+1.d0-gg)/(rk*gg-1.d0/gg)
      func2(gg,rf2de0,psi,rk)=(1.d0+rk*psi)**2*(gg**2-1.d0)-rf2de0
c----------------------------------------------------------------------|

      rk=gm/(gm-1.d0)

      xxmin=1.d0

      do i=i0,i1

        mflg(i)=0

        eede(i)=ee(i)/de(i)
        rf2de(i)=(rx(i)**2)/(de(i)**2)
        psi=func1(xxmin,eede(i),rk)
        fmin=func2(xxmin,rf2de(i),psi,rk)
        psi=func1(xxmax,eede(i),rk)
        fmax=func2(xxmax,rf2de(i),psi,rk)
        if (fmin.lt.-tolf) then
          xx(i)=xxmin
          dxx(i)=xxmax-xxmin
        else if (fmin.gt.tolf) then
          xx(i)=xxmax
          dxx(i)=xxmin-xxmax
        else
          mflg(i)=1
          xx(i)=xxmin
        endif
      enddo

      do mi=1,mix
      do i=i0,i1
      if (mflg(i).eq.0) then

        dxx(i)=dxx(i)*0.5
        xxmid=xx(i)+dxx(i)
        psi=func1(xxmid,eede(i),rk)
        fmid(i)=func2(xxmid,rf2de(i),psi,rk)
        if (fmid(i).le.0.) xx(i)=xxmid
        if (abs(dxx(i)).lt.tolx .or. fmid(i).eq.0.) mflg(i)=1

      endif
      enddo
      enddo

      do i=i0,i1
        gl(i)=xx(i)
        ro(i)=de(i)/gl(i)
        psi=func1(gl(i),eede(i),rk)
        pr(i)=psi*ro(i)
        denom=(1.d0+rk*psi)*gl(i)*de(i)
        vx(i)=rx(i)/denom
      enddo

      return
      end
