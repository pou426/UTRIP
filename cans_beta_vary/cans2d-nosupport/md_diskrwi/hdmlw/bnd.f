c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,ix,jx,roi,pri,vyi)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension roi(ix,jx),pri(ix,jx),vyi(ix,jx)
c----------------------------------------------------------------------|      
      call bdinix(0,margin,ro,roi,ix,jx)
      call bdinix(0,margin,pr,pri,ix,jx)
      call bdcnsx(0,margin,vx,0.d0,ix,jx)
      call bdinix(0,margin,vy,vyi,ix,jx)

      call bdinix(1,margin,ro,roi,ix,jx)
      call bdinix(1,margin,pr,pri,ix,jx)
      call bdcnsx(1,margin,vx,0.d0,ix,jx)
      call bdinix(1,margin,vy,vyi,ix,jx)

      call bdpery(margin,margin,ro,ix,jx)
      call bdpery(margin,margin,pr,ix,jx)
      call bdpery(margin,margin,vx,ix,jx)
      call bdpery(margin,margin,vy,ix,jx)

      return
      end
