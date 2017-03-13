c======================================================================|
      subroutine cooldef(cool,cool0,tecl0,decl0,ro,te,ix,jx)
c======================================================================|
c
c NAME  cooldef
c
c PURPOSE
c    define non-dimensional (approximate) cooling function
c
c OUTPUTS
c    cool(ix,jx) : [double] cooling function
c
c INPUTS
c    cool0: [double] strength of radiative cooling
c    tecl0: [double] temperature for peak of cooling function
c                cooling time is tau_cool=tecl0/cool0
c    decl0: [double] threshold for optically-thick
c    te(ix,jx): [double] temperature
c    ro(ix,jx): [double] density
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension cool(ix,jx)
      dimension ro(ix,jx)
      dimension te(ix,jx)
c----------------------------------------------------------------------|

         a=1.5
         b=2.0
         tmp0=1./(a+b)*log(b/a)
         tmp1=exp(a*tmp0)+exp(-b*tmp0)

      do j=1,jx
      do i=1,ix

         tel0=log10(te(i,j)/tecl0)
         coolfn=0.4*tel0-3.0
     &           +3.0*tmp1/(exp(a*(tel0+tmp0))+exp(-b*(tel0+tmp0)))
         den=decl0*tanh(ro(i,j)/decl0)
         cool(i,j)=cool0*den*10.**coolfn

       enddo
       enddo


      return
      end
