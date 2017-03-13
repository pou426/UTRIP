c======================================================================|
      subroutine bnd(margin,ro,pr,vx,gl,glm,dxm,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix),pr(ix),vx(ix)
      dimension dxm(ix)
c----------------------------------------------------------------------|      
      call bdfrdx(0,margin,ro,dxm,ix)
      call bdfrdx(0,margin,pr,dxm,ix)
      call bdfrdx(0,margin,vx,dxm,ix)
      call bdfrdx(0,margin,gl,dxm,ix)
      call bdfrdx(0,margin-1,glm,dx,ix)

c     call bdfrdx(1,margin,ro,dxm,ix)
c     call bdfrdx(1,margin,pr,dxm,ix)
c     call bdcnsx(1,margin,pr,1.d-3,ix)
c     call bdfrdx(1,margin,vx,dxm,ix)
c     call bdfrdx(1,margin,gl,dxm,ix)
c     call bdfrdx(1,margin,glm,dx,ix)

      call bdfrex(1,margin,ro,ix)
      call bdfrex(1,margin,pr,ix)
      call bdfrex(1,margin,vx,ix)
      call bdfrex(1,margin,gl,ix)
      call bdfrex(1,margin,glm,ix)

      return
      end
