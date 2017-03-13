c======================================================================|
      subroutine bicgstab2(r,rtld,p,v,t,phat,shat,s,work,
     &     rho1,alpha,omega,xx,res,cmat,dmat,margin,ix,jx,mi)
c======================================================================|
c
c NAME  bicgstab2
c
c PURPOSE
c     BiCGStab scheme in order to obtain solution of simultaneous
c     equations, 'Ax=b'.
c    
c INPUTS & OUTPUTS
c     xx(ix,jx): [double] solution, vector "x"
c 
c OUTPUTS
c     res: [double] residue
c     mi: [integer] number of iterations 
c 
c INPUTS
c     ix,jx: [integer] dimension size
c     margin: [integer] size of boundary margins
c     cmat(ix,jx,5): [double] sparse coeffcient matrix, "A"
c     dmat(ix,jx): [double] ILU decomposed diagonal matrix
c     work(ix,jx,7): [double] dimensions for vectors
c     r:    [integer] index for vector "r",  work(i,j,r)
c     rtld: [integer] index for vector "r~", work(i,j,rtrl)
c     p:    [integer] index for vector "p",  work(i,j,p)
c     v:    [integer] index for vector "v",  work(i,j,v)
c     t:    [integer] index for vector "t",  work(i,j,t)
c     phat: [integer] index for vector "p^", work(i,j,p^)
c     shat: [integer] index for vector "s^", work(i,j,s^)
c     s:    [integer] index for vector "s",  work(i,j,s)
c
c HISTORY
c    written 2002-3-23 K. Nakamura
c    vectorized 2006-2-28 K. Nakamura
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,5)
      dimension res(ix,jx)
      dimension work(ix,jx,7),xx(ix,jx),dmat(ix,jx)
      integer   r,rtld,p,v,t,phat,shat,s
      dimension tmpm(ix,jx,5),tmpv(ix),tmpu(ix)
c----------------------------------------------------------------------|

c     rho=(r~,r)
      rho=0.
c      do i=1+margin,ix-margin
c      do j=1+margin,jx-margin
c         rho=rho+work(i,j,rtld)*work(i,j,r)
c      enddo
c      enddo
      do i=1,ix
         tmpv(i)=0.d0
      enddo
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         tmpv(i)=tmpv(i)+work(i,j,rtld)*work(i,j,r)
      enddo
      enddo
      do i=1+margin,ix-margin
         rho=rho+tmpv(i)
      enddo
      
c     p{i}=r+beta*(p{i-1}-omega*v)
      if(mi.gt.1)then
         beta=(rho/rho1)*(alpha/omega)
         do j=margin,jx-margin
         do i=margin,ix-margin
            work(i,j,p)=work(i,j,r)
     &           +beta*(work(i,j,p)-omega*work(i,j,v))
         enddo
         enddo
      else
         do j=1+margin,jx-margin
         do i=1+margin,ix-margin
            work(i,j,p)=work(i,j,r)
         enddo
         enddo
      endif

c     phat = M^{-1} p
      call psolv(work,phat,p,dmat,cmat,margin,ix,jx)

c     v = A phat
c      do i=1+margin,ix-margin
c      do j=1+margin,jx-margin
c         work(i,j,v)=
c     &        +cmat(i,j,2)*work(i-1,j  ,phat)
c     &        +cmat(i,j,4)*work(i  ,j-1,phat)
c     &        +cmat(i,j,1)*work(i  ,j  ,phat)
c     &        +cmat(i,j,5)*work(i  ,j+1,phat)
c     &        +cmat(i,j,3)*work(i+1,j  ,phat)
c      enddo
c      enddo
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
             tmpm(i,j,2)=cmat(i,j,2)*work(i-1,j  ,phat)
             tmpm(i,j,4)=cmat(i,j,4)*work(i  ,j-1,phat)
             tmpm(i,j,1)=cmat(i,j,1)*work(i  ,j  ,phat)
             tmpm(i,j,5)=cmat(i,j,5)*work(i  ,j+1,phat)
             tmpm(i,j,3)=cmat(i,j,3)*work(i+1,j  ,phat)
      enddo
      enddo
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         work(i,j,v)=tmpm(i,j,1)+tmpm(i,j,2)+tmpm(i,j,3)+tmpm(i,j,4)
     &        +tmpm(i,j,5)
      enddo
      enddo

c     alpha
      temp1=0.
c      do i=1+margin,ix-margin
c      do j=1+margin,jx-margin
c         temp1=temp1+work(i,j,rtld)*work(i,j,v)
c      enddo
c      enddo
      do i=1,ix
         tmpv(i)=0.d0
      enddo
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         tmpv(i)=tmpv(i)+work(i,j,rtld)*work(i,j,v)
      enddo
      enddo
      do i=1+margin,ix-margin
         temp1=temp1+tmpv(i)
      enddo
      alpha=rho/temp1

c     s = r - alpha * v
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         work(i,j,s)=work(i,j,r)-alpha*work(i,j,v)
      enddo
      enddo

c     shat = M^{-1} s 
      call psolv(work,shat,s,dmat,cmat,margin,ix,jx)

c     t = A shat
c      do j=1+margin,jx-margin
c      do i=1+margin,ix-margin
c         work(i,j,t)=
c     &        +cmat(i,j,2)*work(i-1,j  ,shat)
c     &        +cmat(i,j,4)*work(i  ,j-1,shat)
c     &        +cmat(i,j,1)*work(i  ,j  ,shat)
c     &        +cmat(i,j,5)*work(i  ,j+1,shat)
c     &        +cmat(i,j,3)*work(i+1,j  ,shat)
c      enddo
c      enddo
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         tmpm(i,j,2)=cmat(i,j,2)*work(i-1,j  ,shat)
         tmpm(i,j,4)=cmat(i,j,4)*work(i  ,j-1,shat)
         tmpm(i,j,1)=cmat(i,j,1)*work(i  ,j  ,shat)
         tmpm(i,j,5)=cmat(i,j,5)*work(i  ,j+1,shat)
         tmpm(i,j,3)=cmat(i,j,3)*work(i+1,j  ,shat)
      enddo
      enddo
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         work(i,j,t)=tmpm(i,j,1)+tmpm(i,j,2)+tmpm(i,j,3)+tmpm(i,j,4)
     &        +tmpm(i,j,5)
      enddo
      enddo

c     omega = (t,s)/(t,t)
      temp2=0.
      temp3=0.
c      do j=1+margin,jx-margin
c      do i=1+margin,ix-margin
c         temp2=temp2+work(i,j,t)*work(i,j,s)
c         temp3=temp3+work(i,j,t)*work(i,j,t)
c      enddo
c      enddo
      do i=1,ix
         tmpv(i)=0.
         tmpu(i)=0.
      enddo
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         tmpv(i)=tmpv(i)+work(i,j,t)*work(i,j,s)
         tmpu(i)=tmpu(i)+work(i,j,t)*work(i,j,t)
      enddo
      enddo
      do i=1+margin,ix-margin
         temp2=temp2+tmpv(i)
         temp3=temp3+tmpu(i)
      enddo

      omega=temp2/temp3

c     x{i} = x{i-1} + alpha*phat + omega*shat 
      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
         xx(i,j)=xx(i,j) + alpha*work(i,j,phat)
     &        + omega*work(i,j,shat)
         work(i,j,r)=work(i,j,s)-omega*work(i,j,t)
      enddo
      enddo
      rho1=rho

c      
      do j=1+margin,jx-1
      do i=1+margin,ix-1
         res(i,j)=work(i,j,s)
      enddo
      enddo
c
      return
      end
