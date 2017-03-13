c======================================================================|
      subroutine p_bbtoaa_c(ay,bz,bx,dzm,dxm,x,margin,ix,jx
     &           ,igx,jgx,ipe,jpe,ipex,jpex)
c======================================================================|
c
c NAME  p_bbtoaa_c
c
c PURPOSE
c    calculate magnetic Vector potential
c        * Cylindrical coordinate, axis-symmetry
c        * Parallel version
c
c OUTPUTS
c    ay(ix,jx): [double] magnetic Vector potential
c
c INPUTS
c    bz(ix,jx),bx(ix,jx): [double] magnetic field
c    dxm(ix),dzm(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dxm(ix),dzm(jx)
      dimension x(ix)
      dimension bx(ix,jx),bz(ix,jx)
      dimension ay(ix,jx)
      dimension ay0(jgx),bx0(jgx),dzm0(jgx),bxg(jgx),dzmg(jgx)
      dimension ayedge(jgx,0:ipex-1),ayedgeg(jgx,0:ipex-1)
      dimension ayoffset(jx)

      include "mpif.h"

c----------------------------------------------------------------------|

      do jg=1,jgx
        bx0(jg)=0.d0
        bxg(jg)=0.d0
        dzm0(jg)=0.d0
        dzmg(jg)=0.d0
      enddo

      if (ipe.eq.0) then
        j0=margin+1
        j1=jx-margin
        if (jpe.eq.0) j0=1
        if (jpe.eq.jpex-1) j1=jx
        do j=j0,j1
           jg=jpe*(jx-2*margin)+j
           bx0(jg)=bx(1,j)
           dzm0(jg)=dzm(j)
        enddo
      endif

      call mpi_allreduce(bx0,bxg,jgx,mpi_double_precision,mpi_sum
     &                      ,mpi_comm_world,merr)
      call mpi_allreduce(dzm0,dzmg,jgx,mpi_double_precision,mpi_sum
     &                      ,mpi_comm_world,merr)

      if (ipe.eq.0) then
        x0=x(margin+1)
        ay0(margin+1)=0.d0
        do jg=margin+2,jgx
          ay0(jg) = ay0(jg-1)
     &          - 0.5*x0*(bxg(jg-1)+bxg(jg  ))*dzmg(jg-1)
        enddo
        do jg=margin,1,-1
          ay0(jg) = ay0(jg+1)
     &          + 0.5*x0*(bxg(jg  )+bxg(jg+1))*dzmg(jg  )
        enddo
      else
        do jg=1,jgx
          ay0(jg) = 0.d0
        enddo
      endif

      do j=1,jx
        jg=jpe*(jx-2*margin)+j
        ay(margin+1,j)=ay0(jg)
      enddo

      do j=1,jx
      do i=margin+2,ix
         ay(i,j) = ay(i-1,j)
     &          + 0.5*(x(i-1)*bz(i-1,j)+x(i)*bz(i  ,j))*dxm(i-1)
      enddo
      do i=margin,1,-1
         ay(i,j) = ay(i+1,j)
     &          - 0.5*(x(i)*bz(i  ,j)+x(i+1)*bz(i+1,j))*dxm(i  )
      enddo
      enddo

      do i=0,ipex-1
      do jg=1,jgx
        ayedge(jg,i)=0.d0
      enddo
      enddo
      do j=margin+1,jx-margin
        jg=jpe*(jx-2*margin)+j
        ayedge(jg,ipe)=ay(ix-margin+1,j)
      enddo
      if (jpe.eq.0) then
       do j=1,margin
         jg=jpe*(jx-2*margin)+j
         ayedge(jg,ipe)=ay(ix-margin+1,j)
       enddo
      endif
      if (jpe.eq.jpex-1) then
       do j=jx-margin+1,jx
         jg=jpe*(jx-2*margin)+j
         ayedge(jg,ipe)=ay(ix-margin+1,j)
       enddo
      endif
      mcnt=jgx*ipex
      call mpi_allreduce(ayedge,ayedgeg,mcnt
     &      ,mpi_double_precision,mpi_sum,mpi_comm_world,merr)

      do j=1,jx
          ayoffset(j)=0.d0
      enddo
      do i=0,ipe-1
      do j=1,jx
        jg=jpe*(jx-2*margin)+j
          ayoffset(j)=ayoffset(j)+ayedgeg(jg,i)
      enddo
      enddo

      do j=1,jx
      do i=1,ix
        ay(i,j)=ay(i,j)+ayoffset(j)
      enddo
      enddo

      do i=1,ix
      do j=1,jx
         ay(i,j) = ay(i,j)/x(i)
      enddo
      enddo


      return
      end
