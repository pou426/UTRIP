c======================================================================|
      subroutine cippre_sh(ro,pr,vx,de,ei,rxm,gm,ix)
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
      dimension ro(ix),pr(ix),vx(ix)
      dimension de(ix),ei(ix),rxm(ix)
c----------------------------------------------------------------------|
         do i=1,ix
           gl=1.d0/sqrt(1.d0-vx(i)**2)
           de(i)=ro(i)*gl
           ei(i)=pr(i)/(gm-1.d0)*gl
         enddo

         do i=1,ix-1
           vxm   =(vx(i)+vx(i+1))/2
           eim   =(ei(i)+ei(i+1))/2
           dem   =(de(i)+de(i+1))/2
           glm   =1.d0/sqrt(1.d0-vxm**2)

           rxm(i)=( eim*gm+dem )*glm*vxm
         enddo

         return
         end
