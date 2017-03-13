c======================================================================|
      subroutine exc_7(margin,da1,da2,da3,da4,da5,da6,da7
     $     ,ix,jx,kx,ipe,jpe,kpe,ipex,jpex,kpex,mper)
c======================================================================|
c   MPI data exchange
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      include 'mpif.h'
      integer mstatus(MPI_STATUS_SIZE)
      dimension da1(ix,jx,kx),da2(ix,jx,kx)
      dimension da3(ix,jx,kx),da4(ix,jx,kx),da5(ix,jx,kx)
      dimension da6(ix,jx,kx),da7(ix,jx,kx)
      parameter(mx=7)
      dimension bufsndx(margin,jx,kx,mx),bufrcvx(margin,jx,kx,mx)
      dimension bufsndy(ix,margin,kx,mx),bufrcvy(ix,margin,kx,mx)
      dimension bufsndz(ix,jx,margin,mx),bufrcvz(ix,jx,margin,mx)
      dimension mper(3)
c----------------------------------------------------------------------|
      kperflag=mper(3)
      jperflag=mper(2)
      iperflag=mper(1)
c----------------------------------------------------------------------|
      mmxx=margin*jx*kx*mx
      mmyx=ix*margin*kx*mx
      mmzx=ix*jx*margin*mx
      mpe=ipex*jpex*kpe+ipex*jpe+ipe
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

      do k=1,kx
      do j=1,jx
      do i=1,margin
         bufsndx(i,j,k,1)=da1(margin+i,j,k)
         bufsndx(i,j,k,2)=da2(margin+i,j,k)
         bufsndx(i,j,k,3)=da3(margin+i,j,k)
         bufsndx(i,j,k,4)=da4(margin+i,j,k)
         bufsndx(i,j,k,5)=da5(margin+i,j,k)
         bufsndx(i,j,k,6)=da6(margin+i,j,k)
         bufsndx(i,j,k,7)=da7(margin+i,j,k)
      enddo
      enddo
      enddo

      call mpi_sendrecv
     &    (bufsndx,mmxx,mpi_double_precision,mped ,0
     &    ,bufrcvx,mmxx,mpi_double_precision,mpeu ,0
     &    ,mpi_comm_world,mstatus,merr)

      if (mpeu.ne.mpi_proc_null) then
      do k=1,kx
      do j=1,jx
      do i=1,margin
         da1(ix-margin+i,j,k)=bufrcvx(i,j,k,1)
         da2(ix-margin+i,j,k)=bufrcvx(i,j,k,2)
         da3(ix-margin+i,j,k)=bufrcvx(i,j,k,3)
         da4(ix-margin+i,j,k)=bufrcvx(i,j,k,4)
         da5(ix-margin+i,j,k)=bufrcvx(i,j,k,5)
         da6(ix-margin+i,j,k)=bufrcvx(i,j,k,6)
         da7(ix-margin+i,j,k)=bufrcvx(i,j,k,7)
      enddo
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

      do k=1,kx
      do j=1,jx
      do i=1,margin
         bufsndx(i,j,k,1)=da1(ix-2*margin+i,j,k)
         bufsndx(i,j,k,2)=da2(ix-2*margin+i,j,k)
         bufsndx(i,j,k,3)=da3(ix-2*margin+i,j,k)
         bufsndx(i,j,k,4)=da4(ix-2*margin+i,j,k)
         bufsndx(i,j,k,5)=da5(ix-2*margin+i,j,k)
         bufsndx(i,j,k,6)=da6(ix-2*margin+i,j,k)
         bufsndx(i,j,k,7)=da7(ix-2*margin+i,j,k)
      enddo
      enddo
      enddo

      call mpi_sendrecv
     &    (bufsndx,mmxx,mpi_double_precision,mpeu ,1
     &    ,bufrcvx,mmxx,mpi_double_precision,mped ,1
     &     ,mpi_comm_world,mstatus,merr)

      if (mped.ne.mpi_proc_null) then
      do k=1,kx
      do j=1,jx
      do i=1,margin
         da1(i,j,k)=bufrcvx(i,j,k,1)
         da2(i,j,k)=bufrcvx(i,j,k,2)
         da3(i,j,k)=bufrcvx(i,j,k,3)
         da4(i,j,k)=bufrcvx(i,j,k,4)
         da5(i,j,k)=bufrcvx(i,j,k,5)
         da6(i,j,k)=bufrcvx(i,j,k,6)
         da7(i,j,k)=bufrcvx(i,j,k,7)
      enddo
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

      do k=1,kx
      do j=1,margin
      do i=1,ix
         bufsndy(i,j,k,1)=da1(i,margin+j,k)
         bufsndy(i,j,k,2)=da2(i,margin+j,k)
         bufsndy(i,j,k,3)=da3(i,margin+j,k)
         bufsndy(i,j,k,4)=da4(i,margin+j,k)
         bufsndy(i,j,k,5)=da5(i,margin+j,k)
         bufsndy(i,j,k,6)=da6(i,margin+j,k)
         bufsndy(i,j,k,7)=da7(i,margin+j,k)
      enddo
      enddo
      enddo

      call mpi_sendrecv
     &    (bufsndy,mmyx,mpi_double_precision,mped ,10
     &    ,bufrcvy,mmyx,mpi_double_precision,mpeu ,10
     &    ,mpi_comm_world,mstatus,merr)

      if (mpeu.ne.mpi_proc_null) then
      do k=1,kx
      do j=1,margin
      do i=1,ix
         da1(i,jx-margin+j,k)=bufrcvy(i,j,k,1)
         da2(i,jx-margin+j,k)=bufrcvy(i,j,k,2)
         da3(i,jx-margin+j,k)=bufrcvy(i,j,k,3)
         da4(i,jx-margin+j,k)=bufrcvy(i,j,k,4)
         da5(i,jx-margin+j,k)=bufrcvy(i,j,k,5)
         da6(i,jx-margin+j,k)=bufrcvy(i,j,k,6)
         da7(i,jx-margin+j,k)=bufrcvy(i,j,k,7)
      enddo
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

      do k=1,kx
      do j=1,margin
      do i=1,ix
         bufsndy(i,j,k,1)=da1(i,jx-2*margin+j,k)
         bufsndy(i,j,k,2)=da2(i,jx-2*margin+j,k)
         bufsndy(i,j,k,3)=da3(i,jx-2*margin+j,k)
         bufsndy(i,j,k,4)=da4(i,jx-2*margin+j,k)
         bufsndy(i,j,k,5)=da5(i,jx-2*margin+j,k)
         bufsndy(i,j,k,6)=da6(i,jx-2*margin+j,k)
         bufsndy(i,j,k,7)=da7(i,jx-2*margin+j,k)
      enddo
      enddo
      enddo

      call mpi_sendrecv
     &    (bufsndy,mmyx,mpi_double_precision,mpeu ,11
     &    ,bufrcvy,mmyx,mpi_double_precision,mped ,11
     &                 ,mpi_comm_world,mstatus,merr)

      if (mped.ne.mpi_proc_null) then
      do k=1,kx
      do j=1,margin
      do i=1,ix
         da1(i,j,k)=bufrcvy(i,j,k,1)
         da2(i,j,k)=bufrcvy(i,j,k,2)
         da3(i,j,k)=bufrcvy(i,j,k,3)
         da4(i,j,k)=bufrcvy(i,j,k,4)
         da5(i,j,k)=bufrcvy(i,j,k,5)
         da6(i,j,k)=bufrcvy(i,j,k,6)
         da7(i,j,k)=bufrcvy(i,j,k,7)
      enddo
      enddo
      enddo
      endif

 2000 continue

c----------------------------------------------------------------------|
c     from PE(myrank) to PE(myrank-1) for new da(ix)
c     z-direction (x-y plane) from mpe to mpe-ipex*jpex
c----------------------------------------------------------------------|
      if (kpex.eq.1) goto 3000

      mpeu = mpe+ipex*jpex
      mped = mpe-ipex*jpex
      if (kperflag.eq.0) then
         if (kpe.eq.kpex-1) mpeu = mpi_proc_null
         if (kpe.eq.0     ) mped = mpi_proc_null
      else
         if (kpe.eq.kpex-1) mpeu = mpe-(kpex-1)*jpex*ipex
         if (kpe.eq.0     ) mped = mpe+(kpex-1)*jpex*ipex
      endif

      do k=1,margin
      do j=1,jx
      do i=1,ix
         bufsndz(i,j,k,1)=da1(i,j,margin+k)
         bufsndz(i,j,k,2)=da2(i,j,margin+k)
         bufsndz(i,j,k,3)=da3(i,j,margin+k)
         bufsndz(i,j,k,4)=da4(i,j,margin+k)
         bufsndz(i,j,k,5)=da5(i,j,margin+k)
         bufsndz(i,j,k,6)=da6(i,j,margin+k)
         bufsndz(i,j,k,7)=da7(i,j,margin+k)
      enddo
      enddo
      enddo

      call mpi_sendrecv
     &    (bufsndz,mmzx,mpi_double_precision,mped ,20
     &    ,bufrcvz,mmzx,mpi_double_precision,mpeu ,20
     &    ,mpi_comm_world,mstatus,merr)

      if (mpeu.ne.mpi_proc_null) then
      do k=1,margin
      do j=1,jx
      do i=1,ix
         da1(i,j,kx-margin+k)=bufrcvz(i,j,k,1)
         da2(i,j,kx-margin+k)=bufrcvz(i,j,k,2)
         da3(i,j,kx-margin+k)=bufrcvz(i,j,k,3)
         da4(i,j,kx-margin+k)=bufrcvz(i,j,k,4)
         da5(i,j,kx-margin+k)=bufrcvz(i,j,k,5)
         da6(i,j,kx-margin+k)=bufrcvz(i,j,k,6)
         da7(i,j,kx-margin+k)=bufrcvz(i,j,k,7)
      enddo
      enddo
      enddo
      endif

c----------------------------------------------------------------------|
c     from PE(myrank) to PE(myrank+1) for new da(1)
c     z-direction (x-y plane) from mpe to mpe+ipex*jpex
c----------------------------------------------------------------------|

      mpeu = mpe+ipex*jpex
      mped = mpe-ipex*jpex
      if (kperflag.eq.0) then
         if (kpe.eq.kpex-1) mpeu = mpi_proc_null
         if (kpe.eq.0     ) mped = mpi_proc_null
      else
         if (kpe.eq.kpex-1) mpeu = mpe-(kpex-1)*jpex*ipex
         if (kpe.eq.0     ) mped = mpe+(kpex-1)*jpex*ipex
      endif

      do k=1,margin
      do j=1,jx
      do i=1,ix
         bufsndz(i,j,k,1)=da1(i,j,kx-2*margin+k)
         bufsndz(i,j,k,2)=da2(i,j,kx-2*margin+k)
         bufsndz(i,j,k,3)=da3(i,j,kx-2*margin+k)
         bufsndz(i,j,k,4)=da4(i,j,kx-2*margin+k)
         bufsndz(i,j,k,5)=da5(i,j,kx-2*margin+k)
         bufsndz(i,j,k,6)=da6(i,j,kx-2*margin+k)
         bufsndz(i,j,k,7)=da7(i,j,kx-2*margin+k)
      enddo
      enddo
      enddo

      call mpi_sendrecv
     &    (bufsndz,mmzx,mpi_double_precision,mpeu ,21
     &    ,bufrcvz,mmzx,mpi_double_precision,mped ,21
     &                 ,mpi_comm_world,mstatus,merr)

      if (mped.ne.mpi_proc_null) then
      do k=1,margin
      do j=1,jx
      do i=1,ix
         da1(i,j,k)=bufrcvz(i,j,k,1)
         da2(i,j,k)=bufrcvz(i,j,k,2)
         da3(i,j,k)=bufrcvz(i,j,k,3)
         da4(i,j,k)=bufrcvz(i,j,k,4)
         da5(i,j,k)=bufrcvz(i,j,k,5)
         da6(i,j,k)=bufrcvz(i,j,k,6)
         da7(i,j,k)=bufrcvz(i,j,k,7)
      enddo
      enddo
      enddo
      endif

 3000 continue

      return
      end
