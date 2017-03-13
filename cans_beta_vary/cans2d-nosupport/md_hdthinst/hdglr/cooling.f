*     ------------------------------------------------   -*- Fortran -*-
*=====================================================================*
*           solving cooling and heating process                       *
*     Koyama & Inutsuka, 2002, ApJ, 564, L97                          *
*=====================================================================*
      subroutine cooling(ro,pr,dt,dtminT,ix,jx)
      implicit none
      integer ix,jx,i,j
      double precision ro(ix,jx),pr(ix,jx),dt,dtminT
      double precision f,df,n1,T1,x,y,coolfn

      coolfn(y,x)= 1.0d0 - y
     &     *(1e7*exp(-1.184e2/(x+1))+.44272*x**.5*exp(-92e-3/x))

      do j=1,jx
      do i=1,ix
         T1 = pr(i,j)/ro(i,j)
         n1 = ro(i,j)
         if(T1.gt.20) goto 1
*        --------------------
         f  = coolfn(n1,T1)
         if(f.eq.0) goto 1
         df = (coolfn(n1,T1*1.001)-f)/(T1*0.001)
         if(f.ne.0)          dtminT = dmin1(dtminT,T1/abs(f))
c         write(*,*)i,j,T1,dt,f
         T1 = T1 + dt * f !/(1-dt*df)
c         if(T1.lt.1e-3)write(*,*)'T1<1e-3'
         pr(i,j) = n1*T1 
 1    enddo
      enddo

      return
      end
