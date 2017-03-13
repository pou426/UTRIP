c======================================================================|
      subroutine bndsc(margin,dsc,dscm,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension dsc(ix),dscm(ix)
c----------------------------------------------------------------------|      
      call bdfrex(0,margin,dsc,ix)
      call bdfrex(0,margin-1,dscm,ix)

      call bdfrex(1,margin,dsc,ix)
      call bdfrex(1,margin,dscm,ix)

      return
      end
