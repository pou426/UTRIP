c======================================================================|
      subroutine intpvl(qql,qqr,qq,cal,car,dxm,dt,ix)
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

      dimension dxm(ix)

      dimension qq(ix)
      dimension qql(ix),qqr(ix)
      dimension cal(ix),car(ix)
      dimension dqq(ix),rqq(ix)
c----------------------------------------------------------------------|


      do i=1,ix-1
          dqq(i)=(qq(i+1)-qq(i))/dxm(i)
      enddo 

      do i=2,ix-1
         if(dqq(i)*dqq(i-1).gt.0.) then
           rqq(i)=2.0*dqq(i)*dqq(i-1)/(dqq(i)+dqq(i-1))
         else
           rqq(i)=0.0
         endif
      enddo 

      do i=2,ix-1
        if(car(i).ge.0.) then
          qqr(i)=qq(i)+0.5*(dxm(i)-car(i)*dt)*rqq(i)
        else
          qqr(i)=qq(i+1)-0.5*(dxm(i)+car(i)*dt)*rqq(i+1)
        endif

        if(cal(i).ge.0.) then
          qql(i)=qq(i)+0.5*(dxm(i)-cal(i)*dt)*rqq(i)
        else
          qql(i)=qq(i+1)-0.5*(dxm(i)+cal(i)*dt)*rqq(i+1)
        endif
      enddo 


      return
      end
