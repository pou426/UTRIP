c======================================================================|
      subroutine bicgstab1(r,rtld,p,v,t,phat,shat,s,work,xx
     &     ,src,margin,ix)
c======================================================================|
c
c NAME  bicgstab1
c
c PURPOSE
c     Initialize the work variables and set useful index to handle the
c     vectors.
c
c OUTPUTS
c     work(ix,7): [double] vectors used in bicgstab2.F
c     xx(ix):     [double] initial guess of vector x (Ax=b)
c     r:    [integer] index for vector "r",  work(ix,r)
c     rtld: [integer] index for vector "r~", work(ix,rtrl)
c     p:    [integer] index for vector "p",  work(ix,p)
c     v:    [integer] index for vector "v",  work(ix,v)
c     t:    [integer] index for vector "t",  work(ix,t)
c     phat: [integer] index for vector "p^", work(ix,p^)
c     shat: [integer] index for vector "s^", work(ix,s^)
c     s:    [integer] index for vector "s",  work(ix,s)
c
c INPUTS
c     ix:      [integer] dimension size
c     margin:  [integer] size of boundary margins
c     src(ix): [double] source vector of heat conduction eq
c
c HISTORY
c    written 2002-3-23 K. Nakamura
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension work(ix,1:7)
      dimension xx(ix)
      dimension src(ix)
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

      do i=1,ix
         xx(i)       =0.
      enddo

c     r=y
      do i=1,margin
        work(i,r)=0.
      enddo
      do i=1+margin,ix-margin
         work(i,r)=src(i)
      enddo

c     r~=r
      do i=1+margin,ix-margin
         work(i,rtld)=work(i,r)
      enddo

      return
      end
