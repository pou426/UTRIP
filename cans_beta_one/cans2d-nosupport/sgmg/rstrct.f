c     -------------------------------------------------------------
      subroutine rstrct(uc,uf,nc)
      integer nc
      double precision uc(-1:nc,-1:nc), uf(-1:2*nc,-1:2*nc)
      integer ic, if, jc, jf
c     --- boundary condition: uf
c     --- begin: periodic boundary condition
      do if=-1,2*nc
         uf(if,-1)=uf(if,2*nc-1)
      enddo
      do jf=-1,2*nc
         uf(-1,jf)=uf(2*nc-1,jf)
      enddo
c     --- end: periodic boundary condition
c
      do jc=0,nc-1
         jf=2*jc
         do ic=0,nc-1
            if=2*ic
            uc(ic,jc)=0.5d0*uf(if,jf)+0.125d0*(uf(if+1,jf)+
     $           uf(if-1,jf)+uf(if,jf+1)+uf(if,jf-1))
         enddo
      enddo
c     --- boundary condition: uc
c     --- begin: periodic boundary condition
      do ic=0,nc
         uc(ic,nc)=uc(ic,0)
      enddo
      do jc=0,nc
         uc(nc,jc)=uc(0,jc)
      enddo
c     --- end: periodic boundary condition
      return
      end
