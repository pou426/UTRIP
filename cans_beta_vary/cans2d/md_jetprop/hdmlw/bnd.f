c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vz,gm,x,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx)
      dimension pr(ix,jx)
      dimension vx(ix,jx)
      dimension vz(ix,jx)
      dimension x(ix)
c----------------------------------------------------------------------|      
c  jet flow

      ro1=0.1d0
      pr1=1.d0/gm
      vx1=0
      vz1=19.d0

c----------------------------------------------------------------------|
      call bdsppx(0,margin,ro,ix,jx)
      call bdsppx(0,margin,pr,ix,jx)
      call bdspnx(0,margin,vx,ix,jx)
      call bdsppx(0,margin,vz,ix,jx)

      call bdfrex(1,margin,ro,ix,jx)
      call bdfrex(1,margin,pr,ix,jx)
      call bdfrex(1,margin,vx,ix,jx)
      call bdfrex(1,margin,vz,ix,jx)

      call bdsppy(0,margin,ro,ix,jx)
      call bdsppy(0,margin,pr,ix,jx)
      call bdsppy(0,margin,vx,ix,jx)
      call bdspny(0,margin,vz,ix,jx)
      do j=1,margin
      do i=1,ix
        if (x(i).le.1.d0) then
          ro(i,j)=ro1
          pr(i,j)=pr1
          vx(i,j)=vx1
          vz(i,j)=vz1
        endif
      enddo
      enddo

      call bdfrey(1,margin,ro,ix,jx)
      call bdfrey(1,margin,pr,ix,jx)
      call bdfrey(1,margin,vx,ix,jx)
      call bdfrey(1,margin,vz,ix,jx)

      return
      end

