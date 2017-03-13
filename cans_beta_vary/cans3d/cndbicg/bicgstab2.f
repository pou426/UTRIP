c======================================================================|
      subroutine bicgstab2(r,rtld,p,v,t,phat,shat,s,work,
     &     rho1,alpha,omega,xx,res,cmat,dmat,margin,ix,jx,kx,mi)
c======================================================================|
c
c NAME  bicgstab2
c
c PURPOSE
c     BiCGStab scheme in order to obtain solution of simultaneous
c     equations, 'Ax=b'.
c    
c INPUTS & OUTPUTS
c     xx(ix,jx,kx): [double] solution, vector "x"
c 
c OUTPUTS
c     res: [double] residue
c     mi: [integer] number of iterations 
c 
c INPUTS
c     ix,jx,kx: [integer] dimension size
c     margin: [integer] size of boundary margins
c     cmat(ix,jx,kx,7): [double] sparse coeffcient matrix, "A"
c     dmat(ix,jx,kx): [double] ILU decomposed diagonal matrix
c     work(ix,jx,kx,7): [double] dimensions for vectors
c     r:    [integer] index for vector "r",  work(i,j,k,r)
c     rtld: [integer] index for vector "r~", work(i,j,k,rtrl)
c     p:    [integer] index for vector "p",  work(i,j,k,p)
c     v:    [integer] index for vector "v",  work(i,j,k,v)
c     t:    [integer] index for vector "t",  work(i,j,k,t)
c     phat: [integer] index for vector "p^", work(i,j,k,p^)
c     shat: [integer] index for vector "s^", work(i,j,k,s^)
c     s:    [integer] index for vector "s",  work(i,j,k,s)
c
c HISTORY
c    written 2004-3-26 K. Nakamura
c    vectorized 2006-2-28 K. Nakamura
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,kx,7)
      dimension res(ix,jx,kx)
      dimension work(ix,jx,kx,7),xx(ix,jx,kx),dmat(ix,jx,kx)
      integer   r,rtld,p,v,t,phat,shat,s
      dimension tmpm(ix,jx,kx,7),tmpv(ix),tmpu(ix)
c----------------------------------------------------------------------|

c     rho=(r~,r)
      rho=0.
c      do k=1+margin,kx-margin
c      do j=1+margin,jx-margin
c      do i=1+margin,ix-margin
c         rho=rho+work(i,j,k,rtld)*work(i,j,k,r)
c      enddo
c      enddo
c      enddo
      do i=1+margin,ix-margin
         tmpv(i)=0.
      enddo
      do k=1+margin,kx-margin
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         tmpv(i)=tmpv(i)+work(i,j,k,rtld)*work(i,j,k,r)
      enddo
      enddo
      enddo
      do i=1+margin,ix-margin
         rho=rho+tmpv(i)
      enddo
      
c     p{i}=r+beta*(p{i-1}-omega*v)
      if(mi.gt.1)then
         beta=(rho/rho1)*(alpha/omega)
         do k=margin,kx-margin
         do j=margin,jx-margin
         do i=margin,ix-margin
            work(i,j,k,p)=work(i,j,k,r)
     &           +beta*(work(i,j,k,p)-omega*work(i,j,k,v))
         enddo
         enddo
         enddo
      else
         do k=1+margin,kx-margin
         do j=1+margin,jx-margin
         do i=1+margin,ix-margin
            work(i,j,k,p)=work(i,j,k,r)
         enddo
         enddo
         enddo
      endif

c     phat = M^{-1} p
      call psolv(work,phat,p,dmat,cmat,margin,ix,jx,kx)

c     v = A phat
c      do k=1+margin,kx-margin
c      do j=1+margin,jx-margin
c      do i=1+margin,ix-margin
c         work(i,j,k,v)=
c     &        +cmat(i,j,k,2)*work(i-1,j  ,k  ,phat)
c     &        +cmat(i,j,k,4)*work(i  ,j-1,k  ,phat)
c     &        +cmat(i,j,k,6)*work(i  ,j  ,k-1,phat)
c     &        +cmat(i,j,k,1)*work(i  ,j  ,k  ,phat)
c     &        +cmat(i,j,k,7)*work(i  ,j  ,k+1,phat)
c     &        +cmat(i,j,k,5)*work(i  ,j+1,k  ,phat)
c     &        +cmat(i,j,k,3)*work(i+1,j  ,k  ,phat)
c      enddo
c      enddo
c      enddo
      do k=1+margin,kx-margin
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         tmpm(i,j,k,2)=cmat(i,j,k,2)*work(i-1,j  ,k  ,phat)
         tmpm(i,j,k,4)=cmat(i,j,k,4)*work(i  ,j-1,k  ,phat)
         tmpm(i,j,k,6)=cmat(i,j,k,6)*work(i  ,j  ,k-1,phat)
         tmpm(i,j,k,1)=cmat(i,j,k,1)*work(i  ,j  ,k  ,phat)
         tmpm(i,j,k,7)=cmat(i,j,k,7)*work(i  ,j  ,k+1,phat)
         tmpm(i,j,k,5)=cmat(i,j,k,5)*work(i  ,j+1,k  ,phat)
         tmpm(i,j,k,3)=cmat(i,j,k,3)*work(i+1,j  ,k  ,phat)
      enddo
      enddo
      enddo
      do k=1+margin,kx-margin
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         work(i,j,k,v)=tmpm(i,j,k,1)+tmpm(i,j,k,2)+tmpm(i,j,k,3)
     &        +tmpm(i,j,k,4)+tmpm(i,j,k,5)+tmpm(i,j,k,6)+tmpm(i,j,k,7)
      enddo
      enddo
      enddo

c     alpha
      temp1=0.
c      do k=1+margin,kx-margin
c      do j=1+margin,jx-margin
c      do i=1+margin,ix-margin
c         temp1=temp1+work(i,j,k,rtld)*work(i,j,k,v)
c      enddo
c      enddo
c      enddo
      do i=1+margin,ix-margin
         tmpv(i)=0.
      enddo
      do k=1+margin,kx-margin
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         tmpv(i)=tmpv(i)+work(i,j,k,rtld)*work(i,j,k,v)
      enddo
      enddo
      enddo
      do i=1+margin,ix-margin
         temp1=temp1+tmpv(i)
      enddo
      
      alpha=rho/temp1

c     s = r - alpha * v
      do k=1+margin,kx-margin
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         work(i,j,k,s)=work(i,j,k,r)-alpha*work(i,j,k,v)
      enddo
      enddo
      enddo

c     shat = M^{-1} s 
      call psolv(work,shat,s,dmat,cmat,margin,ix,jx,kx)

c     t = A shat
c      do k=1+margin,kx-margin
c      do j=1+margin,jx-margin
c      do i=1+margin,ix-margin
c         work(i,j,k,t)=
c     &        +cmat(i,j,k,2)*work(i-1,j  ,k  ,shat)
c     &        +cmat(i,j,k,4)*work(i  ,j-1,k  ,shat)
c     &        +cmat(i,j,k,6)*work(i  ,j  ,k-1,shat)
c     &        +cmat(i,j,k,1)*work(i  ,j  ,k  ,shat)
c     &        +cmat(i,j,k,7)*work(i  ,j  ,k+1,shat)
c     &        +cmat(i,j,k,5)*work(i  ,j+1,k  ,shat)
c     &        +cmat(i,j,k,3)*work(i+1,j  ,k  ,shat)
c      enddo
c      enddo
c      enddo
      do k=1+margin,kx-margin
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         tmpm(i,j,k,2) = cmat(i,j,k,2)*work(i-1,j  ,k  ,shat)
         tmpm(i,j,k,4) = cmat(i,j,k,4)*work(i  ,j-1,k  ,shat)
         tmpm(i,j,k,6) = cmat(i,j,k,6)*work(i  ,j  ,k-1,shat)
         tmpm(i,j,k,1) = cmat(i,j,k,1)*work(i  ,j  ,k  ,shat)
         tmpm(i,j,k,7) = cmat(i,j,k,7)*work(i  ,j  ,k+1,shat)
         tmpm(i,j,k,5) = cmat(i,j,k,5)*work(i  ,j+1,k  ,shat)
         tmpm(i,j,k,3) = cmat(i,j,k,3)*work(i+1,j  ,k  ,shat)
      enddo
      enddo
      enddo
      do k=1+margin,kx-margin
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         work(i,j,k,t) = tmpm(i,j,k,1)+tmpm(i,j,k,2)+tmpm(i,j,k,3)
     &        +tmpm(i,j,k,4)+tmpm(i,j,k,5)+tmpm(i,j,k,6)+tmpm(i,j,k,7)
      enddo
      enddo
      enddo

c     omega = (t,s)/(t,t)
      temp2=0.
      temp3=0.
c      do k=1+margin,kx-margin
c      do j=1+margin,jx-margin
c      do i=1+margin,ix-margin
c         temp2=temp2+work(i,j,k,t)*work(i,j,k,s)
c         temp3=temp3+work(i,j,k,t)*work(i,j,k,t)
c      enddo
c      enddo
c      enddo
      do i=1+margin,ix-margin
         tmpv(i)=0.
         tmpu(i)=0.
      enddo
      do k=1+margin,kx-margin
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         tmpv(i)=tmpv(i)+work(i,j,k,t)*work(i,j,k,s)
         tmpu(i)=tmpu(i)+work(i,j,k,t)*work(i,j,k,t)
      enddo
      enddo
      enddo
      do i=1+margin,ix-margin
         temp2=temp2+tmpv(i)
         temp3=temp3+tmpu(i)
      enddo
      omega=temp2/temp3

c     x{i} = x{i-1} + alpha*phat + omega*shat 
      do k=1+margin,kx-margin
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         xx(i,j,k)=xx(i,j,k) + alpha*work(i,j,k,phat)
     &        + omega*work(i,j,k,shat)
         work(i,j,k,r)=work(i,j,k,s)-omega*work(i,j,k,t)
      enddo
      enddo
      enddo
      rho1=rho

c      
      do k=1+margin,kx-1
      do j=1+margin,jx-1
      do i=1+margin,ix-1
         res(i,j,k)=work(i,j,k,s)
      enddo
      enddo
      enddo
c
      return
      end
