c======================================================================|
      subroutine cippost_sm(de,ei,rxm,ry,ro,pr,vx,vy,gm,ix)
c======================================================================|
c
c NAME  prtote
c
c PURPOSE
c    calculation of temperature
c
c OUTPUTS
c    te(ix): [double] temperature
c
c INPUTS
c    ro(ix): [double] density
c    pr(ix): [double] pressure
c    gm: [double] polytropic index gamma
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension ro(ix),pr(ix),vx(ix),vy(ix)
      dimension de(ix),ei(ix),rxm(ix),ry(ix)
c----------------------------------------------------------------------|

      do i=2,ix
        rx=(rxm(i-1)+rxm(i))/2
        rt = sqrt( rx**2 + ry(i)**2 + ( de(i)+gm*ei(i) )**2 )
        vx(i) = rx/rt
        vy(i) = ry(i)/rt
        gl = 1./sqrt(1.-vx(i)**2-vy(i)**2)
        pr(i) = ei(i)/gl*(gm-1.d0)
        ro(i) = de(i)/gl
      enddo

      return
      end
