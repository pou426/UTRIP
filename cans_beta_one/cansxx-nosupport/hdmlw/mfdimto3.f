c======================================================================|
      subroutine mfdimto3(mfdim,mfdim3)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension mfdim(3),mfdim3(3,3)

      mfdim3(1,1)=mfdim(1)
      mfdim3(2,1)=mfdim(2)
      mfdim3(3,1)=mfdim(3)

      mfdim3(1,2)=mfdim(2)
      mfdim3(2,2)=mfdim(3)
      mfdim3(3,2)=mfdim(1)

      mfdim3(1,3)=mfdim(3)
      mfdim3(2,3)=mfdim(1)
      mfdim3(3,3)=mfdim(2)

      return
      end
