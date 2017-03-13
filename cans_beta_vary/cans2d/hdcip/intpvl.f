c======================================================================|
      subroutine intpvl(qql,qqr,qq,cal,car,dxm,dym,dt,ix,jx,mdir)
c======================================================================|
c
c NAME  intpvl
c
c PURPOSE
c    van Leer interpolation
c
c OUTPUTS
c    qql(ix), qqr(ix): [double] interpolated results.
c             defined at the left (right) side of grid boundary
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    qq(ix): [double] data of variable to be interpolated
c    cal(ix), car(ix) : [double] characteristic speed 
c                of leftward (rightward) propagating wave
c    dxm(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dxm(ix),dym(jx)

      dimension qq(ix,jx)
      dimension qql(ix,jx),qqr(ix,jx)
      dimension cal(ix,jx),car(ix,jx)
      dimension dqq(ix,jx),rqq(ix,jx)
c----------------------------------------------------------------------|
      if (mdir.eq.1) then
        ip=1
        jp=0
      else
        ip=0
        jp=1
      endif


      do j=1,jx-jp
      do i=1,ix-ip
          ds=dxm(i)*ip+dym(j)*jp
          dqq(i,j)=(qq(i+ip,j+jp)-qq(i,j))/ds
      enddo 
      enddo 

      do j=1+jp,jx-jp
      do i=1+ip,ix-ip
         if(dqq(i,j)*dqq(i-ip,j-jp).gt.0.) then
           rqq(i,j)=2.0*dqq(i,j)*dqq(i-ip,j-jp)
     &              /(dqq(i,j)+dqq(i-ip,j-jp))
         else
           rqq(i,j)=0.0
         endif
      enddo 
      enddo 

      do j=1+jp,jx-jp
      do i=1+ip,ix-ip
        ds=dxm(i)*ip+dym(j)*jp
        if(car(i,j).ge.0.) then
          qqr(i,j)=qq(i,j)+0.5*(ds-car(i,j)*dt)*rqq(i,j)
        else
          qqr(i,j)=qq(i+ip,j+jp)-0.5*(ds+car(i,j)*dt)*rqq(i+ip,j+jp)
        endif

        if(cal(i,j).ge.0.) then
          qql(i,j)=qq(i,j)+0.5*(ds-cal(i,j)*dt)*rqq(i,j)
        else
          qql(i,j)=qq(i+ip,j+jp)-0.5*(ds+cal(i,j)*dt)*rqq(i+ip,j+jp)
        endif
      enddo 
      enddo 


      return
      end
