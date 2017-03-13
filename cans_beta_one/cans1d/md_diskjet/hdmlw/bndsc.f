c======================================================================|
      subroutine bndsc(margin,sc,dsc,scm,dscm,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension sc(ix),dsc(ix)
      dimension scm(ix),dscm(ix)
c----------------------------------------------------------------------|      
c     call bdfrex(0,margin,sc,ix)
      call bdfrex(0,margin,dsc,ix)
c     call bdfrex(0,margin-1,scm,ix)
      call bdfrex(0,margin-1,dscm,ix)

c     call bdfrex(1,margin,sc,ix)
      call bdfrex(1,margin,dsc,ix)
c     call bdfrex(1,margin,scm,ix)
      call bdfrex(1,margin,dscm,ix)

      call bdsppx(0,margin,sc,ix)
      call bdspnx(0,margin,dsc,ix)
      call bdsmpx(0,margin-1,scm,ix)
      call bdsmnx(0,margin-1,dscm,ix)

      call bdsppx(1,margin,sc,ix)
      call bdspnx(1,margin,dsc,ix)
      call bdsmpx(1,margin,scm,ix)
      call bdsmnx(1,margin,dscm,ix)


      return
      end
