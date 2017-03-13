c======================================================================|
      subroutine intpvl(qql,qqr,qq,cal,car,dxm,dym,dzm,dt
     &                     ,ix,jx,kx,mdir)
c======================================================================|
c
c NAME  intpvl
c
c PURPOSE
c    van Leer interpolation
c
c OUTPUTS
c    qql(ix,jx,kx), qqr(ix,jx,kx): [double] interpolated results.
c             defined at the left (right) side of grid boundary
c
c INPUTS
c    NOTE: ??m(ix,jx,kx) is the variable array defined at grid bounds
c
c    qq(ix,jx,kx): [double] data of variable to be interpolated
c    cal(ix,jx,kx), car(ix,jx,kx) : [double] characteristic speed 
c                of leftward (rightward) propagating wave
c    dxm(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix,jx,kx: [integer] dimension size
c
c HISTORY
c    written 2003-6-1 K. Takahashi based on T. Yokoyama's code
c
c----------------------------------------------------------------------|

      implicit real*8 (a-h,o-z)

      dimension dxm(ix),dym(jx),dzm(kx)

      dimension qq(ix,jx,kx)
      dimension qql(ix,jx,kx),qqr(ix,jx,kx)
      dimension cal(ix,jx,kx),car(ix,jx,kx)
      dimension dqq(ix,jx,kx),rqq(ix,jx,kx)
c----------------------------------------------------------------------|
      if (mdir.eq.1) then
        ip=1
        jp=0
        kp=0
      endif
      if (mdir.eq.2) then
        ip=0
        jp=1
        kp=0
      endif
      if (mdir.eq.3) then
        ip=0
        jp=0
        kp=1
      endif


      do k=1,kx-kp
      do j=1,jx-jp
      do i=1,ix-ip
          ds=dxm(i)*ip+dym(j)*jp+dzm(k)*kp
          dqq(i,j,k)=(qq(i+ip,j+jp,k+kp)-qq(i,j,k))/ds
      enddo
      enddo
      enddo

      do k=1+kp,kx-kp
      do j=1+jp,jx-jp
      do i=1+ip,ix-ip
         if(dqq(i,j,k)*dqq(i-ip,j-jp,k-kp).gt.0.) then
           rqq(i,j,k)=2.0*dqq(i,j,k)*dqq(i-ip,j-jp,k-kp)
     &                  /(dqq(i,j,k)+dqq(i-ip,j-jp,k-kp))
         else
           rqq(i,j,k)=0.0
         endif
      enddo
      enddo
      enddo

      do k=1+kp,kx-kp
      do j=1+jp,jx-jp
      do i=1+ip,ix-ip
        ds=dxm(i)*ip+dym(j)*jp+dzm(k)*kp
        if(car(i,j,k).ge.0.) then
          qqr(i,j,k)=qq(i,j,k)
     &              +0.5*(ds-car(i,j,k)*dt)*rqq(i,j,k)
        else
          qqr(i,j,k)=qq(i+ip,j+jp,k+kp)
     &              -0.5*(ds+car(i,j,k)*dt)*rqq(i+ip,j+jp,k+kp)
        endif

        if(cal(i,j,k).ge.0.) then
          qql(i,j,k)=qq(i,j,k)
     &              +0.5*(ds-cal(i,j,k)*dt)*rqq(i,j,k)
        else
          qql(i,j,k)=qq(i+ip,j+jp,k+kp)
     &              -0.5*(ds+cal(i,j,k)*dt)*rqq(i+ip,j+jp,k+kp)
        endif
      enddo
      enddo
      enddo


      return
      end
