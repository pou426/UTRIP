c======================================================================|
      subroutine mlwsf(u2,dt,s,ux0,ux1,uy0,uy1,uz0,uz1,ix,jx,kx,mfdim)
c======================================================================|
c
c NAME  mlwsrcf
c
c PURPOSE
c    apply source term for second half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    u2(ix,jx,kx): [double] variation in this step
c
c OUTPUTS
c    None
c
c INPUTS
c    s(ix,jx,kx) : [double] source term
c    ux0(ix) : [double] 0.5*dx(i)/dxm(i)
c    ux1(ix) : [double] 0.5*dx(i-1)/dxm(i)
c    dt: [double] delta time 
c    ix,jx,kx: [integer] dimension size
c 
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension mfdim(3)
      dimension u2(ix,jx,kx),s(ix,jx,kx)
      dimension ux0(ix),ux1(ix)
      dimension uy0(jx),uy1(jx)
      dimension uz0(kx),uz1(kx)
c----------------------------------------------------------------------|      

c     do k=2,kx-1
c     do j=2,jx-1
c     do i=2,ix-1
c       sh = uz0(k)*(uy0(j)*(ux0(i)*s(i-1,j-1,k-1)+ux1(i)*s(i,j-1,k-1))
c    &             + uy1(j)*(ux0(i)*s(i-1,j,k-1)  +ux1(i)*s(i,j,k-1)  ))
c    &     + uz1(k)*(uy0(j)*(ux0(i)*s(i-1,j-1,k)  +ux1(i)*s(i,j-1,k)  )
c    &             + uy1(j)*(ux0(i)*s(i-1,j,k)    +ux1(i)*s(i,j,k)    ))
c       u2(i,j,k)= u2(i,j,k)+0.5d0*dt*sh
c     enddo
c     enddo
c     enddo

c----------------------------------------------------------------------|      
c  3D-xyz

      if (mfdim(1).eq.1.and.mfdim(2).eq.1.and.mfdim(3).eq.1) then
      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
        u2(i,j,k) =u2(i,j,k) +0.5d0*dt
     &    *(uz0(k)*(uy0(j)*(ux0(i)*s(i-1,j-1,k-1)+ux1(i)*s(i,j-1,k-1))
     &             +uy1(j)*(ux0(i)*s(i-1,j  ,k-1)+ux1(i)*s(i,j  ,k-1)))
     &     +uz1(k)*(uy0(j)*(ux0(i)*s(i-1,j-1,k)  +ux1(i)*s(i,j-1,k)  )
     &             +uy1(j)*(ux0(i)*s(i-1,j  ,k)  +ux1(i)*s(i,j  ,k)  )))
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|      
c  2D-xy

      else if (mfdim(1).eq.1.and.mfdim(2).eq.1.and.mfdim(3).eq.0) then
      do k=1,kx
      do j=2,jx-1
      do i=2,ix-1
        u2(i,j,k) =u2(i,j,k) +0.5d0*dt
     &            *(uy0(j)*(ux0(i)*s(i-1,j-1,k  )+ux1(i)*s(i,j-1,k  ))
     &             +uy1(j)*(ux0(i)*s(i-1,j  ,k  )+ux1(i)*s(i,j  ,k  )))
      enddo
      enddo
      enddo

c  2D-zx

      else if (mfdim(1).eq.1.and.mfdim(2).eq.0.and.mfdim(3).eq.1) then
      do k=2,kx-1
      do j=1,jx
      do i=2,ix-1
        u2(i,j,k) =u2(i,j,k) +0.5d0*dt
     &    *(uz0(k)*((ux0(i)*s(i-1,j  ,k-1)+ux1(i)*s(i,j  ,k-1)))
     &     +uz1(k)*((ux0(i)*s(i-1,j  ,k)  +ux1(i)*s(i,j  ,k)  )))
      enddo
      enddo
      enddo

c  2D-yz

      else if (mfdim(1).eq.0.and.mfdim(2).eq.1.and.mfdim(3).eq.1) then
      do k=2,kx-1
      do j=2,jx-1
      do i=1,ix
        u2(i,j,k) =u2(i,j,k) +0.5d0*dt
     &    *(uz0(k)*(uy0(j)*s(i,j-1,k-1)
     &             +uy1(j)*s(i,j  ,k-1))
     &     +uz1(k)*(uy0(j)*s(i,j-1,k)  
     &             +uy1(j)*s(i,j  ,k)  ))
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|      
c  1D-x

      else if (mfdim(1).eq.1.and.mfdim(2).eq.0.and.mfdim(3).eq.0) then
      do k=1,kx
      do j=1,jx
      do i=2,ix-1
        u2(i,j,k) =u2(i,j,k) +0.5d0*dt
     &    * (ux0(i)*s(i-1,j  ,k)  +ux1(i)*s(i,j  ,k)  )
      enddo
      enddo
      enddo

c  1D-y

      else if (mfdim(1).eq.0.and.mfdim(2).eq.1.and.mfdim(3).eq.0) then
      do k=1,kx
      do j=2,jx-1
      do i=1,ix
        u2(i,j,k) =u2(i,j,k) +0.5d0*dt
     &    *(uy0(j)*s(i,j-1,k)  
     &     +uy1(j)*s(i,j  ,k)  )
      enddo
      enddo
      enddo

c  1D-z

      else if (mfdim(1).eq.0.and.mfdim(2).eq.0.and.mfdim(3).eq.1) then
      do k=2,kx-1
      do j=1,jx
      do i=1,ix
        u2(i,j,k) =u2(i,j,k) +0.5d0*dt
     &    *(uz0(k)*s(i,j  ,k-1)
     &     +uz1(k)*s(i,j  ,k)  )
      enddo
      enddo
      enddo

      endif
      
      return
      end
