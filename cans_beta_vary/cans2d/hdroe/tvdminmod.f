      subroutine tvdminmod(mdir,da,daw,ix,jx)
      implicit double precision (a-h,o-z)

      dimension da(ix,jx)
      dimension daw(ix,jx,2)

c----------------------------------------------------------------------|
c     define limiter functions
c     1. minmod limiter

      flmt(a,b)=max(0.0d0,min(b*sign(1.0d0,a),abs(a)))*sign(1.0d0,a)
c----------------------------------------------------------------------|

      if (mdir.eq.1) then
        do j=2,jx-2
        do i=2,ix-2
           daw(i,j,1)=da(i,j)  
     &               +0.5*flmt(da(i+1,j)-da(i,j),da(i,j)-da(i-1,j))
           daw(i,j,2)=da(i+1,j)
     &               -0.5*flmt(da(i+1,j)-da(i,j),da(i+2,j)-da(i+1,j))
        enddo
        enddo
      endif

      if (mdir.eq.2) then
        do j=2,jx-2
        do i=2,ix-2
           daw(i,j,1)=da(i,j)  
     &               +0.5*flmt(da(i,j+1)-da(i,j),da(i,j)-da(i,j-1))
           daw(i,j,2)=da(i,j+1)
     &               -0.5*flmt(da(i,j+1)-da(i,j),da(i,j+2)-da(i,j+1))
        enddo
        enddo
      endif

      return
      end
