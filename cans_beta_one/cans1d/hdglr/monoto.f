      subroutine monoto(dx,var,grad,ix)

      implicit double precision (a-h,o-z)

      dimension dx(ix),var(ix),grad(ix)
      dimension dx1inv(ix),dx2inv(ix)

      do i=1,ix-1
         dx1inv(i) = 0.5d0/dx(i)
         dx2inv(i) = 2.0d0/dx(i)
      enddo

      do i=2,ix-1
         agrad = (var(i+1) - var(i-1))*dx1inv(i)
         xgrad = (var(i+1) - var(i))  *dx2inv(i)
         wgrad = (var(i)   - var(i-1))*dx2inv(i)
         if( xgrad*wgrad .gt. 0.0d0 ) then
            grad(i) = wgrad*xgrad/(wgrad + xgrad)
         else
            grad(i) = 0.0d0
         endif
      enddo

      return
      end

