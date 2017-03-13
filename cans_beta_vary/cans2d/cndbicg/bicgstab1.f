c======================================================================|
      subroutine bicgstab1(r,rtld,p,v,t,phat,shat,s,work,xx,src
     &                     ,margin,ix,jx)
c======================================================================|
c
c NAME  bicgstab1
c
c PURPOSE
c     Initialize the work variables and set useful index to handle the
c     vectors.
c
c OUTPUTS
c     work(ix,jx,7): [double] vectors used in bicgstab2.F
c     xx(ix,jx):     [double] initial guess of vector x (Ax=b)
c     r:    [integer] index for vector "r",  work(ix,jx,r)
c     rtld: [integer] index for vector "r~", work(ix,jx,rtrl)
c     p:    [integer] index for vector "p",  work(ix,jx,p)
c     v:    [integer] index for vector "v",  work(ix,jx,v)
c     t:    [integer] index for vector "t",  work(ix,jx,t)
c     phat: [integer] index for vector "p^", work(ix,jx,p^)
c     shat: [integer] index for vector "s^", work(ix,jx,s^)
c     s:    [integer] index for vector "s",  work(ix,jx,s)
c
c INPUTS
c     ix,jx:   [integer] dimension size
c     margin:  [integer] size of boundary margins
c     src(ix,jx): [double] source vector of heat conduction eq
c
c HISTORY
c    written 2002-3-23 K. Nakamura
c    vectorized 2006-2-28 K. Nakamura
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension work(ix,jx,7)
      dimension xx(ix,jx)
      dimension src(ix,jx)
      integer r,rtld,p,v,t,phat,shat,s
c----------------------------------------------------------------------|
c     A x = y (cmat x = src)
      r    = 1
      rtld = 2
      p    = 3
      v    = 4
      t    = 5
      phat = 6
      shat = 7
      s    = 1

      do j=1,jx
      do i=1,ix
         xx(i,j)       = 0.
      enddo
      enddo

      do j=1+margin,jx-margin
      do i=1+margin,ix-margin
c     r=y
         work(i,j,r)   = src(i,j)

c     r~=r
         work(i,j,rtld)= src(i,j)
      enddo
      enddo
c
      return
      end
