c======================================================================|
      subroutine p_bbtoaa(az,bx,by,dxm,dym,margin,ix,jx
     &            ,igx,jgx,ipe,jpe,ipex,jpex)
c======================================================================|
c
c NAME  p_bbtoaa
c
c PURPOSE
c    calculate magnetic Vector potential
c        * Parallel version
c
c OUTPUTS
c    az(ix,jx): [double] magnetic Vector potential
c
c INPUTS
c    bx(ix,jx),by(ix,jx): [double] magnetic field
c    dxm(ix),dym(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dxm(ix),dym(jx)
      dimension bx(ix,jx),by(ix,jx)
      dimension az(ix,jx)
      dimension az0(jgx),bx0(jgx),dym0(jgx),bxg(jgx),dymg(jgx)
      dimension azedge(jgx,0:ipex-1),azedgeg(jgx,0:ipex-1)
      dimension azoffset(jx)

      include "mpif.h"

c----------------------------------------------------------------------|

      do jg=1,jgx
        bx0(jg)=0.d0
        bxg(jg)=0.d0
        dym0(jg)=0.d0
        dymg(jg)=0.d0
      enddo

      if (ipe.eq.0) then
        j0=margin+1
        j1=jx-margin
        if (jpe.eq.0) j0=1
        if (jpe.eq.jpex-1) j1=jx
        do j=j0,j1
           jg=jpe*(jx-2*margin)+j
           bx0(jg)=bx(margin+1,j)
           dym0(jg)=dym(j)
        enddo
      endif

      call mpi_allreduce(bx0,bxg,jgx,mpi_double_precision,mpi_sum
     &                      ,mpi_comm_world,merr)
      call mpi_allreduce(dym0,dymg,jgx,mpi_double_precision,mpi_sum
     &                      ,mpi_comm_world,merr)

      if (ipe.eq.0) then
        az0(margin+1)=0.d0
        do jg=margin+2,jgx
          az0(jg) = az0(jg-1)
     &          + 0.5*(bxg(jg-1)+bxg(jg  ))*dymg(jg-1)
        enddo
        do jg=margin,1,-1
          az0(jg) = az0(jg+1)
     &          - 0.5*(bxg(jg  )+bxg(jg+1))*dymg(jg  )
        enddo
      else
        do jg=1,jgx
          az0(jg) = 0.d0
        enddo
      endif

      do j=1,jx
        jg=jpe*(jx-2*margin)+j
        az(margin+1,j)=az0(jg)
      enddo

      do j=1,jx
      do i=margin+2,ix
         az(i,j) = az(i-1,j)
     &          - 0.5*(by(i-1,j)+by(i  ,j))*dxm(i-1)
      enddo
      do i=margin,1,-1
         az(i,j) = az(i+1,j)
     &          + 0.5*(by(i  ,j)+by(i+1,j))*dxm(i  )
      enddo
      enddo

      do i=0,ipex-1
      do jg=1,jgx
        azedge(jg,i)=0.d0
      enddo
      enddo
      do j=margin+1,jx-margin
        jg=jpe*(jx-2*margin)+j
        azedge(jg,ipe)=az(ix-margin+1,j)
      enddo
      if (jpe.eq.0) then
       do j=1,margin
         jg=jpe*(jx-2*margin)+j
         azedge(jg,ipe)=az(ix-margin+1,j)
       enddo
      endif
      if (jpe.eq.jpex-1) then
       do j=jx-margin+1,jx
         jg=jpe*(jx-2*margin)+j
         azedge(jg,ipe)=az(ix-margin+1,j)
       enddo
      endif
      mcnt=jgx*ipex
      call mpi_allreduce(azedge,azedgeg,mcnt
     &      ,mpi_double_precision,mpi_sum,mpi_comm_world,merr)

      do j=1,jx
          azoffset(j)=0.d0
      enddo

      do i=0,ipe-1
      do j=1,jx
        jg=jpe*(jx-2*margin)+j
          azoffset(j)=azoffset(j)+azedgeg(jg,i)
      enddo
      enddo

      do j=1,jx
      do i=1,ix
        az(i,j)=az(i,j)+azoffset(j)
      enddo
      enddo

      return
      end
