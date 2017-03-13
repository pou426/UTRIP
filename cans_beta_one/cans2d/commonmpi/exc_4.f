c======================================================================|
      subroutine exc_4(margin,da1,da2,da3,da4
     $     ,ix,jx,ipe,jpe,ipex,jpex,mper)
c======================================================================|
c   MPI data exchange
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      include 'mpif.h'
      integer mstatus(MPI_STATUS_SIZE)
      dimension da1(ix,jx),da2(ix,jx)
      dimension da3(ix,jx),da4(ix,jx)
      parameter(mx=4)
      dimension bufsndx(margin,jx,mx),bufrcvx(margin,jx,mx)
      dimension bufsndy(ix,margin,mx),bufrcvy(ix,margin,mx)
      dimension mper(2)
c----------------------------------------------------------------------|
      jperflag=mper(2)
      iperflag=mper(1)
c----------------------------------------------------------------------|
      mmxx=margin*jx*mx
      mmyx=ix*margin*mx
      mpe=ipex*jpe+ipe
c----------------------------------------------------------------------|
c     from PE(myrank) to PE(myrank-1) for new da(ix)
c     x-direction (y-z plane) from mpe to mpe-1
c----------------------------------------------------------------------|
      if (ipex.eq.1) goto 1000

      mpeu = mpe+1
      mped = mpe-1
      if (iperflag.eq.0) then
         if (ipe.eq.ipex-1) mpeu = mpi_proc_null
         if (ipe.eq.0     ) mped = mpi_proc_null
      else
         if (ipe.eq.ipex-1) mpeu = mpe-(ipex-1)
         if (ipe.eq.0     ) mped = mpe+(ipex-1)
      endif

      do j=1,jx
      do i=1,margin
         bufsndx(i,j,1)=da1(margin+i,j)
         bufsndx(i,j,2)=da2(margin+i,j)
         bufsndx(i,j,3)=da3(margin+i,j)
         bufsndx(i,j,4)=da4(margin+i,j)
      enddo
      enddo

      call mpi_sendrecv
     &    (bufsndx,mmxx,mpi_double_precision,mped ,0
     &    ,bufrcvx,mmxx,mpi_double_precision,mpeu ,0
     &    ,mpi_comm_world,mstatus,merr)

      if (mpeu.ne.mpi_proc_null) then
      do j=1,jx
      do i=1,margin
         da1(ix-margin+i,j)=bufrcvx(i,j,1)
         da2(ix-margin+i,j)=bufrcvx(i,j,2)
         da3(ix-margin+i,j)=bufrcvx(i,j,3)
         da4(ix-margin+i,j)=bufrcvx(i,j,4)
      enddo
      enddo
      endif

c----------------------------------------------------------------------|
c     from PE(myrank) to PE(myrank+1) for new da(1)
c     x-direction (y-z plane) from mpe to mpe+1
c----------------------------------------------------------------------|

      mpeu = mpe+1
      mped = mpe-1
      if (iperflag.eq.0) then
         if (ipe.eq.ipex-1) mpeu = mpi_proc_null
         if (ipe.eq.0     ) mped = mpi_proc_null
      else
         if (ipe.eq.ipex-1) mpeu = mpe-(ipex-1)
         if (ipe.eq.0     ) mped = mpe+(ipex-1)
      endif

      do j=1,jx
      do i=1,margin
         bufsndx(i,j,1)=da1(ix-2*margin+i,j)
         bufsndx(i,j,2)=da2(ix-2*margin+i,j)
         bufsndx(i,j,3)=da3(ix-2*margin+i,j)
         bufsndx(i,j,4)=da4(ix-2*margin+i,j)
      enddo
      enddo

      call mpi_sendrecv
     &    (bufsndx,mmxx,mpi_double_precision,mpeu ,1
     &    ,bufrcvx,mmxx,mpi_double_precision,mped ,1
     &     ,mpi_comm_world,mstatus,merr)

      if (mped.ne.mpi_proc_null) then
      do j=1,jx
      do i=1,margin
         da1(i,j)=bufrcvx(i,j,1)
         da2(i,j)=bufrcvx(i,j,2)
         da3(i,j)=bufrcvx(i,j,3)
         da4(i,j)=bufrcvx(i,j,4)
      enddo
      enddo
      endif

 1000 continue

c----------------------------------------------------------------------|
c     from PE(myrank) to PE(myrank-1) for new da(ix)
c     y-direction (x-z plane) from mpe to mpe-ipex
c----------------------------------------------------------------------|
      if (jpex.eq.1) goto 2000

      mpeu = mpe+ipex
      mped = mpe-ipex
      if (jperflag.eq.0) then
         if (jpe.eq.jpex-1) mpeu = mpi_proc_null
         if (jpe.eq.0     ) mped = mpi_proc_null
      else
         if (jpe.eq.jpex-1) mpeu = mpe-(jpex-1)*ipex
         if (jpe.eq.0     ) mped = mpe+(jpex-1)*ipex
      endif

      do j=1,margin
      do i=1,ix
         bufsndy(i,j,1)=da1(i,margin+j)
         bufsndy(i,j,2)=da2(i,margin+j)
         bufsndy(i,j,3)=da3(i,margin+j)
         bufsndy(i,j,4)=da4(i,margin+j)
      enddo
      enddo

      call mpi_sendrecv
     &    (bufsndy,mmyx,mpi_double_precision,mped ,10
     &    ,bufrcvy,mmyx,mpi_double_precision,mpeu ,10
     &    ,mpi_comm_world,mstatus,merr)

      if (mpeu.ne.mpi_proc_null) then
      do j=1,margin
      do i=1,ix
         da1(i,jx-margin+j)=bufrcvy(i,j,1)
         da2(i,jx-margin+j)=bufrcvy(i,j,2)
         da3(i,jx-margin+j)=bufrcvy(i,j,3)
         da4(i,jx-margin+j)=bufrcvy(i,j,4)
      enddo
      enddo
      endif

c----------------------------------------------------------------------|
c     from PE(myrank) to PE(myrank+1) for new da(1)
c     y-direction (x-z plane) from mpe to mpe+ipex
c----------------------------------------------------------------------|

      mpeu = mpe+ipex
      mped = mpe-ipex
      if (jperflag.eq.0) then
         if (jpe.eq.jpex-1) mpeu = mpi_proc_null
         if (jpe.eq.0     ) mped = mpi_proc_null
      else
         if (jpe.eq.jpex-1) mpeu = mpe-(jpex-1)*ipex
         if (jpe.eq.0     ) mped = mpe+(jpex-1)*ipex
      endif

      do j=1,margin
      do i=1,ix
         bufsndy(i,j,1)=da1(i,jx-2*margin+j)
         bufsndy(i,j,2)=da2(i,jx-2*margin+j)
         bufsndy(i,j,3)=da3(i,jx-2*margin+j)
         bufsndy(i,j,4)=da4(i,jx-2*margin+j)
      enddo
      enddo

      call mpi_sendrecv
     &    (bufsndy,mmyx,mpi_double_precision,mpeu ,11
     &    ,bufrcvy,mmyx,mpi_double_precision,mped ,11
     &                 ,mpi_comm_world,mstatus,merr)

      if (mped.ne.mpi_proc_null) then
      do j=1,margin
      do i=1,ix
         da1(i,j)=bufrcvy(i,j,1)
         da2(i,j)=bufrcvy(i,j,2)
         da3(i,j)=bufrcvy(i,j,3)
         da4(i,j)=bufrcvy(i,j,4)
      enddo
      enddo
      endif

 2000 continue


      return
      end
