c======================================================================|
      subroutine rtbis_sh(vx,gm,de,ee,rx,ix,i0,i1)
c======================================================================|
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension vx(ix)
      dimension de(ix),ee(ix),rx(ix)
      dimension fmid(ix),rtbis(ix),dxx(ix),mflg(ix)
c----------------------------------------------------------------------|
      func(vx0,de0,ee0,rx0)=rx0/de0*(1.-(gm-1.)/gm*(1.-vx0**2))
     &        -(ee0/de0+1.-(gm-1.)/gm*sqrt(1.-vx0**2))*vx0
c----------------------------------------------------------------------|

      acc=1.d-10
      mx=1000

      vx1=-0.99999999999999999999999999999999999999999
      vx2=+0.99999999999999999999999999999999999999999

      do i=i0,i1
        mflg(i)=0
        fmid(i)=func(vx2,de(i),ee(i),rx(i))
        f=func(vx1,de(i),ee(i),rx(i))
        if (f.lt.0.) then
          rtbis(i)=vx1
          dxx(i)=vx2-vx1
        else
          rtbis(i)=vx2
          dxx(i)=vx1-vx2
        endif
      enddo

      do m=1,mx
      do i=i0,i1
      if (mflg(i).eq.0) then

        dxx(i)=dxx(i)*0.5
        xmid=rtbis(i)+dxx(i)
        fmid(i)=func(xmid,de(i),ee(i),rx(i))
        if (fmid(i).le.0.) rtbis(i)=xmid
        if (abs(dxx(i)).lt.acc .or. fmid(i).eq.0.) mflg(i)=1

      endif
      enddo
      enddo

      do i=i0,i1
        vx(i)=rtbis(i)
      enddo

      return
      end
