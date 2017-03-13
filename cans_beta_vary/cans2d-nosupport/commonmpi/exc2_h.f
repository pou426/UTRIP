c======================================================================|
      subroutine exc2_h(margin,mx,ro,pr,vx,vy,ix,jx,myrank,npe
     & ,mpi_win,bufrcv)
c======================================================================|
c   MPI data exchange
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      include 'mpif.h'
      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx)
      dimension bufsnd1(ix,margin,mx),bufsnd2(ix,margin,mx)
      dimension bufrcv(ix,margin,mx,2)
c----------------------------------------------------------------------|
c  from PE(myrank) to PE(myrank-1) for new da(ix)
c----------------------------------------------------------------------|

      mright=  myrank+1
      mleft =  myrank-1
      if (myrank.eq.npe-1) mright  = mpi_proc_null
      if (myrank.eq.0    ) mleft   = mpi_proc_null

      do j=1,margin
      do i=1,ix
        bufsnd1(i,j,1)=ro(i,margin+j)
        bufsnd1(i,j,2)=pr(i,margin+j)
        bufsnd1(i,j,3)=vx(i,margin+j)
        bufsnd1(i,j,4)=vy(i,margin+j)
      enddo
      enddo

      call mpi_put(bufsnd1,ix*margin*mx,mpi_double_precision,mleft,
     &                 ,0
     &                 ,ix*margin*mx,mpi_double_precision
     &                 ,mpi_win)

c----------------------------------------------------------------------|
c  from PE(myrank) to PE(myrank+1) for new da(1)
c----------------------------------------------------------------------|

      mright=  myrank+1
      mleft =  myrank-1
      if (myrank.eq.npe-1) mright  = mpi_proc_null
      if (myrank.eq.0    ) mleft   = mpi_proc_null

      do j=1,margin
      do i=1,ix
        bufsnd2(i,j,1)=ro(i,jx-2*margin+j)
        bufsnd2(i,j,2)=pr(i,jx-2*margin+j)
        bufsnd2(i,j,3)=vx(i,jx-2*margin+j)
        bufsnd2(i,j,4)=vy(i,jx-2*margin+j)
      enddo
      enddo

      call mpi_put(bufsnd2,ix*margin*mx,mpi_double_precision,mright,
     &                 ,ix*margin*mx*8
     &                 ,ix*margin*mx,mpi_double_precision
     &                 ,mpi_win)

c----------------------------------------------------------------------|
c  recieve
c----------------------------------------------------------------------|

      call mpi_win_fence(0,mpi_win)

      if (myrank.ne.npe-1) then
      do j=1,margin
      do i=1,ix
        ro(i,jx-margin+j)=bufrcv(i,j,1,1)
        pr(i,jx-margin+j)=bufrcv(i,j,2,1)
        vx(i,jx-margin+j)=bufrcv(i,j,3,1)
        vy(i,jx-margin+j)=bufrcv(i,j,4,1)
      enddo
      enddo
      endif

      if (myrank.ne.0) then
      do j=1,margin
      do i=1,ix
        ro(i,j)=bufrcv(i,j,1,2)
        pr(i,j)=bufrcv(i,j,2,2)
        vx(i,j)=bufrcv(i,j,3,2)
        vy(i,j)=bufrcv(i,j,4,2)
      enddo
      enddo
      endif


      return
      end
