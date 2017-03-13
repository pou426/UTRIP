c======================================================================|
      subroutine cipbnd(margin,te,vxm,vy,by,rodx,tedx,vxdxm,vyrdx
     &  ,ro,dx,dxm,rr,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension te(ix),vxm(ix),vy(ix),by(ix),ro(ix)
      dimension dx(ix),dxm(ix),rr(ix)
      dimension rodx(ix),tedx(ix),vxdxm(ix),vyrdx(ix)
c----------------------------------------------------------------------|      
      call bdspnx(0,margin,rodx,ix)
      call bdsppx(0,margin,te,ix)
      call bdspnx(0,margin,tedx,ix)
      call bdsmnx(0,margin-1,vxm,ix)
      call bdsmpx(0,margin-1,vxdxm,ix)
      call bdsppx(0,margin,vy,ix)
      call bdspnx(0,margin,vyrdx,ix)
      call bdsppx(0,margin,by,ix)


      zero=0.0d0
      call bdfrex(1,margin,te,ix)
      call bdfrex(1,margin-1,vxm,ix)
      call bdfrex(1,margin,vy,ix)
      call bdfrex(1,margin,by,ix)

c     call bdcnsx(1,margin-1,rodx,zero,ix)
c     call bdcnsx(1,margin-1,tedx,zero,ix)
c     call bdcnsx(1,margin-1,vxdxm,zero,ix)
c     call bdcnsx(1,margin-1,vyrdx,zero,ix)

      do i0=1,margin-1
        i=ix-margin+i0
        rodx(i)=(ro(i+1)-ro(i-1))/dx(i)/2
        tedx(i)=(te(i+1)-te(i-1))/dx(i)/2
        vyrdx(i)=(vy(i+1)*rr(i+1)-vy(i-1)*rr(i-1))/dx(i)/2
        vxdxm(i)=(vxm(i+1)-vxm(i-1))/dxm(i)/2
      enddo

c     call bdfrex(1,margin,rodx,ix)
c     call bdfrex(1,margin,tedx,ix)
c     call bdfrex(1,margin,vxdxm,ix)
c     call bdfrex(1,margin,vyrdx,ix)


      return
      end
