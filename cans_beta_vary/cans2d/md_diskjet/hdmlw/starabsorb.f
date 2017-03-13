c======================================================================|
      subroutine starabsorb(ro,pr,vx,vy,vz,bx,by,bz,ay,x,z,ix,jx)
c======================================================================|
c     apply boundary condition
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx),vz(ix,jx)
      dimension bx(ix,jx),by(ix,jx),bz(ix,jx)
      dimension ay(ix,jx)
      dimension x(ix),z(jx)
c----------------------------------------------------------------------|

      rstar=0.2d0
      wdamp=0.004
      rofloor=1.d-4
      prfloor=1.d-4

      do j=1,jx
      do i=1,ix
         rr=sqrt(x(i)**2+z(j)**2)
         if (rr.le.rstar) then
            da=0.1*(1-tanh((rr-rstar+2.5*wdamp)/wdamp))
            ro(i,j)=ro(i,j)-da*(ro(i,j)-rofloor)
            pr(i,j)=pr(i,j)-da*(pr(i,j)-prfloor)
            vx(i,j)=vx(i,j)-da*vx(i,j)
            vy(i,j)=vy(i,j)-da*vy(i,j)
            vz(i,j)=vz(i,j)-da*vz(i,j)
            bx(i,j)=bx(i,j)-da*bx(i,j)
            by(i,j)=by(i,j)-da*by(i,j)
            bz(i,j)=bz(i,j)-da*bz(i,j)
            ay(i,j)=ay(i,j)-da*ay(i,j)
         endif
      enddo
      enddo

      return
      end
