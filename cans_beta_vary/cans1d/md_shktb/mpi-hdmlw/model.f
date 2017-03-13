c======================================================================|
      subroutine model(ro,pr,vx,gm,margin,x,dxm,ix,mf_params,igx,ipe)
c======================================================================|
      implicit real*8 (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension ro(ix),pr(ix),vx(ix)

      dimension xg(igx),dxmg(igx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=1.4d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
c      dxm,x

      dx0=1.d0/real(igx-margin*2)
      do i=1,igx
         dxmg(i)=dx0
      enddo
       
      izero=igx/2
      xg(izero)=-dxmg(izero)/2
      do i=izero+1,igx
         xg(i) = xg(i-1)+dxmg(i-1)
      enddo
      do i=izero-1,1,-1
         xg(i) = xg(i+1)-dxmg(i)
      enddo

      do i=1,ix
         ig=ipe*(ix-2*margin)+i
         x(i)=xg(ig)
         dxm(i)=dxmg(ig)
      enddo
c----------------------------------------------------------------------|
c     initial condition
c----------------------------------------------------------------------|
      ro1=0.125d0
      pr1=0.1d0
      vx0=0.d0
      vx1=0.d0

      do i=1,ix
         if (x(i).le.0.d0) then
           ro(i)  = 1.d0
           pr(i)  = 1.d0
           vx(i)  = vx0
         else
           ro(i)  = ro1
           pr(i)  = pr1
           vx(i)  = vx1
         endif
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'ro1',ro1)
      call dacputparamd(mf_params,'pr1',pr1)
      call dacputparamd(mf_params,'vx0',vx0)
      call dacputparamd(mf_params,'vx1',vx1)

      
      return
      end
