c======================================================================|
      subroutine cipbnd(margin,de,ei,rxm,dedx,eidx,rxdxm,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension de(ix),ei(ix),rxm(ix)
      dimension dedx(ix),eidx(ix),rxdxm(ix)
c----------------------------------------------------------------------|      
      call bdfrex(0,margin,de,ix)
      call bdfrex(0,margin,dedx,ix)
      call bdfrex(0,margin,ei,ix)
      call bdfrex(0,margin,eidx,ix)
      call bdfrex(0,margin-1,rxm,ix)
      call bdfrex(0,margin-1,rxdxm,ix)

      call bdfrex(1,margin,de,ix)
      call bdfrex(1,margin,dedx,ix)
      call bdfrex(1,margin,ei,ix)
      call bdfrex(1,margin,eidx,ix)
      call bdfrex(1,margin,rxm,ix)
      call bdfrex(1,margin,rxdxm,ix)

      return
      end
