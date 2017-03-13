*=====================================================================*
*           solving thermal conduction                                *
*     Koyama & Inutsuka, 2002, ApJ, 564, L97                          *
*=====================================================================*
      subroutine CONDUCT(ro,pr,dt,dtminC,gm,ix,jx,dx,dy)
      implicit none 
      integer ix,jx,i,j
      double precision ro(ix,jx),pr(ix,jx),dtminC,gm,dx(ix),dy(jx)
      double precision cs(ix,jx),fluxx(ix,jx),fluxy(ix,jx)
      INTEGER dfdx,dfdy
      double precision  dt,dxdxi,critc,Kc
      parameter (Kc = 0.000625d0)

      do j=1,jx
      do i=1,ix
         cs(i,j)=sqrt(gm*pr(i,j)/ro(i,j))
      enddo
      enddo

      do j = 1, jx
      do i = 2, ix
         fluxx(i,j) = (cs(i,j) + cs(i-1,j))*(pr(i,j) - pr(i-1,j))
     &                     /(dx(i) + dx(i-1))
      enddo
      enddo

      do j = 2, jx
      do i = 1, ix
         fluxy(i,j) = (cs(i,j) + cs(i,j-1))*(pr(i,j) - pr(i,j-1))
     &                     /(dy(j) + dy(j-1))
      enddo
      enddo

      do j = 2, jx-1
      do i = 2, ix-1
c         write(*,*)'p=',i,j,pr(i,j)
         dfdx = (fluxx(i+1,j)-fluxx(i,j))/dx(i)
         dfdy = (fluxy(i,j+1)-fluxy(i,j))/dy(j)
         pr(i,j)=pr(i,j)+dt*(gm-1)*Kc*(dfdx+dfdy)/ro(i,j)
c         if(pr(i,j)/ro(i,j).lt.1e-3)then
c            write(*,*)i,j,pr(i,j),dt,(gm-1)*Kc*(dfdx+dfdy)/ro(i,j)
c            stop
c            endif
         dtminC=dmin1(dtminC,dx(i)**2*ro(i,j)/((gm-1)*Kc*cs(i,j)))
         dtminC=dmin1(dtminC,dy(j)**2*ro(i,j)/((gm-1)*Kc*cs(i,j)))
      enddo
      enddo
c      write(*,*)dt,dtminC
      return
      end



