c======================================================================|
      subroutine bicgstab2(r,rtld,p,v,t,phat,shat,s,work,
     &     rho1,alpha,omega,xx,res,cmat,dmat,margin,ix,mi)
c======================================================================|
c
c NAME  bicgstab2
c
c PURPOSE
C     BiCGStab scheme in order to obtain solution of simultaneous
C     equations, 'Ax=b'.
c
c INPUTS & OUTPUTS
c     xx(ix): [double] solution, vector "x"
c
c OUTPUTS
c     res: [double] residue
c     mi: [integer] number of iterations
c
c INPUTS
c     ix: [integer] dimension size
c     margin: [integer] size of boundary margins
c     cmat(ix,3): [double] sparse coeffcient matrix, "A"
c     dmat(ix): [double] ILU decomposed diagonal matrix
c     work(ix,7): [double] dimensions for vectors
c     r:    [integer] index for vector "r",  work(i,r)
c     rtld: [integer] index for vector "r~", work(i,rtrl)
c     p:    [integer] index for vector "p",  work(i,p)
c     v:    [integer] index for vector "v",  work(i,v)
c     t:    [integer] index for vector "t",  work(i,t)
c     phat: [integer] index for vector "p^", work(i,p^)
c     shat: [integer] index for vector "s^", work(i,s^)
c     s:    [integer] index for vector "s",  work(i,s)
c
c HISTORY
c    written 2002-3-23 K. Nakamura
c
c----------------------------------------------------------------------|
      
      implicit double precision (a-h,o-z)
      dimension cmat(ix,1:3)
      dimension res(ix)
      dimension work(ix,1:7),xx(ix),dmat(ix)
      integer r,rtld,p,v,t,phat,shat,s

C----------------------------------------------------------------------|

c     rho=(r~,r)
      rho=0.
      do i=1+margin,ix-margin
         rho=rho+work(i,rtld)*work(i,r)
      enddo

c     p{i}=r+beta*(p{i-1}-omega*v)
      if(mi.gt.1)then
         beta=(rho/rho1)*(alpha/omega)
         do i=1+margin,ix-margin
            work(i,p)=work(i,r)+beta*(work(i,p)-omega*work(i,v))
         enddo
      else
         do i=1+margin,ix-margin
            work(i,p)=work(i,r)
         enddo
      endif

c     phat = M^{-1} p
      call psolv(work,phat,p,dmat,cmat,margin,ix)

c     v = A phat
      do i=1+margin,ix-margin
         work(i,v)=
     &        +cmat(i,2)*work(i-1,phat)
     &        +cmat(i,1)*work(i  ,phat)
     &        +cmat(i,3)*work(i+1,phat)
      enddo

c     alpha
      temp1=0.
      do i=1+margin,ix-margin
         temp1=temp1+work(i,rtld)*work(i,v)
      enddo
      alpha=rho/temp1

c     s = r - alpha * v
      do i=1+margin,ix-margin
         work(i,s)=work(i,r)-alpha*work(i,v)
      enddo

c     shat = M^{-1} s 
      call psolv(work,shat,s,dmat,cmat,margin,ix)

c     t = A shat
      do i=1+margin,ix-margin
         work(i,t)=
     &        +cmat(i,2)*work(i-1,shat)
     &        +cmat(i,1)*work(i  ,shat)
     &        +cmat(i,3)*work(i+1,shat)
      enddo

c     omega = (t,s)/(t,t)
      temp2=0.
      temp3=0.
      do i=1+margin,ix-margin
         temp2=temp2+work(i,t)*work(i,s)
         temp3=temp3+work(i,t)*work(i,t)
      enddo
      omega=temp2/temp3

c     x{i} = x{i-1} + alpha*phat + omega*shat 
      do i=1+margin,ix-margin
         xx(i)=xx(i) + alpha*work(i,phat) + omega*work(i,shat)
         work(i,r)=work(i,s)-omega*work(i,t)
      enddo

      rho1=rho
      
      do i=1+margin,ix-margin
         res(i)=work(i,s)
      enddo

      return
      end
