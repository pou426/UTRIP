c======================================================================|
      subroutine exc_3(margin,da1,da2,da3
     $     ,ix,ipe,ipex,mper)
c======================================================================|
c   MPI data exchange
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      include 'mpif.h'
      integer mstatus(MPI_STATUS_SIZE)
      dimension da1(ix),da2(ix)
      dimension da3(ix)
      parameter(mx=3)
      dimension bufsndx(margin,mx),bufrcvx(margin,mx)
c----------------------------------------------------------------------|
      iperflag=mper
c----------------------------------------------------------------------|
      mmxx=margin*mx
      mpe=ipe
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

      do i=1,margin
         bufsndx(i,1)=da1(margin+i)
         bufsndx(i,2)=da2(margin+i)
         bufsndx(i,3)=da3(margin+i)
      enddo

      call mpi_sendrecv
     &    (bufsndx,mmxx,mpi_double_precision,mped ,0
     &    ,bufrcvx,mmxx,mpi_double_precision,mpeu ,0
     &    ,mpi_comm_world,mstatus,merr)

      if (mpeu.ne.mpi_proc_null) then
      do i=1,margin
         da1(ix-margin+i)=bufrcvx(i,1)
         da2(ix-margin+i)=bufrcvx(i,2)
         da3(ix-margin+i)=bufrcvx(i,3)
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

      do i=1,margin
         bufsndx(i,1)=da1(ix-2*margin+i)
         bufsndx(i,2)=da2(ix-2*margin+i)
         bufsndx(i,3)=da3(ix-2*margin+i)
      enddo

      call mpi_sendrecv
     &    (bufsndx,mmxx,mpi_double_precision,mpeu ,1
     &    ,bufrcvx,mmxx,mpi_double_precision,mped ,1
     &     ,mpi_comm_world,mstatus,merr)

      if (mped.ne.mpi_proc_null) then
      do i=1,margin
         da1(i)=bufrcvx(i,1)
         da2(i)=bufrcvx(i,2)
         da3(i)=bufrcvx(i,3)
      enddo
      endif

 1000 continue


      return
      end
