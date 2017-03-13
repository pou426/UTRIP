      subroutine monoto(dx,var,grad,ix,jx,id,jd)

      implicit double precision (a-h,o-z)

      dimension dx(ix,jx),var(ix,jx),grad(ix,jx)
      dimension dx1inv(ix,jx),dx2inv(ix,jx)

      do j=1,jx-jd
      do i=1,ix-id
         dx1inv(i,j) = 0.5d0/dx(i,j)
         dx2inv(i,j) = 2.0d0/dx(i,j)
      enddo
      enddo

      do j=1+jd,jx-jd
      do i=1+id,ix-id
         agrad = (var(i+id,j+jd) - var(i-id,j-jd))*dx1inv(i,j)
         xgrad = (var(i+id,j+jd) - var(i,j))  *dx2inv(i,j)
         wgrad = (var(i,j)   - var(i-id,j-jd))*dx2inv(i,j)
         if( xgrad*wgrad .gt. 0.0d0 ) then
            grad(i,j) = wgrad*xgrad/(wgrad + xgrad)
         else
            grad(i,j) = 0.0d0
         endif
      enddo
      enddo

      return
      end
