c     ------------------------------------------------------------------
      subroutine addint(uf,uc,res,nf)
      integer nf
      double precision res(-1:nf,-1:nf),uc(-1:nf/2,-1:nf/2),
     1                 uf(-1:nf,-1:nf)
c     uses interp
      integer i,j
      call interp(res,uc,nf)
      do j=0,nf
         do i=0,nf
            uf(i,j)=uf(i,j)+res(i,j)
         enddo
      enddo
      return
      end
