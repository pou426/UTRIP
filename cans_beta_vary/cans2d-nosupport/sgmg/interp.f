c     ----------------------------------------------------------------
      subroutine interp(uf,uc,nf)
      integer nf
      double precision uc(-1:nf/2,-1:nf/2),uf(-1:nf,-1:nf)
      integer ic,if,jc,jf,nc
      nc=nf/2
      do jf=0,nf,2
         jc=jf/2
         do if=0,nf,2
            ic=if/2
            uf(if,jf)=uc(ic,jc)
         enddo
      enddo
      do jf=0,nf,2
         do if=1,nf-1,2
            uf(if,jf)=0.5d0*(uf(if+1,jf)+uf(if-1,jf))
         enddo
      enddo
      do jf=1,nf-1,2
         do if=0,nf
            uf(if,jf)=0.5d0*(uf(if,jf+1)+uf(if,jf-1))
         enddo
      enddo
      return
      end
