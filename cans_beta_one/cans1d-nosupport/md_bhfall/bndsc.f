c======================================================================|
      subroutine bndsc(margin,hh,hhm,gg,ggm,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension hh(ix),gg(ix)
      dimension hhm(ix),ggm(ix)
c----------------------------------------------------------------------|      
c     call bdsppx(0,margin,hh,ix)
c     call bdspnx(0,margin,gg,ix)
c     call bdsmpx(0,margin-1,hhm,ix)
c     call bdsmnx(0,margin-1,ggm,ix)

c     call bdsppx(1,margin,hh,ix)
c     call bdspnx(1,margin,gg,ix)
c     call bdsmpx(1,margin,hhm,ix)
c     call bdsmnx(1,margin,ggm,ix)

      return
      end

