c======================================================================|
      subroutine bnd(margin,pr,ix,jx
     &           ,ipe,jpe,ipex,jpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension pr(ix,jx)
c----------------------------------------------------------------------|      
      if (ipe.eq.0) then
      call bdsppx(0,margin,pr,ix,jx)
      endif

      if (ipe.eq.ipex-1) then
      call bdsppx(1,margin,pr,ix,jx)
      endif

      if (jpe.eq.0) then
      call bdsppy(0,margin,pr,ix,jx)
      endif

      if (jpe.eq.jpex-1) then
      call bdsppy(1,margin,pr,ix,jx)
      endif


      return
      end
