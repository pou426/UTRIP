c=====================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
      parameter (ipex=1,jpex=2,mpex=ipex*jpex)
      parameter (margin=1)
      parameter (ix=32/ipex+2*margin)
      parameter (jx=32/jpex+2*margin)

      dimension x(ix),xm(ix),dx(ix),dxm(ix)
      dimension y(jx),ym(jx),dy(jx),dym(jx)

      dimension ro(ix,jx),pr(ix,jx)

      include "mpif.h"
      character cno*4
      dimension mper(2)

913   format (1x,' write    ','step=',i8,' time=',e10.3,' nd =',i3)
915   format (1x,' stop     ','step=',i8,' time=',e10.3)
925   format (1x,' ns= ',i10,' mi= ',i5,' err= ',e14.6)
c======================================================================|
c     prologue
c======================================================================|
      merr  = 0
      mcont=0
c----------------------------------------------------------------------|
c     set parameters controling finalization and data-output
      tend=1.0d0
      dtout=0.1d0
      nstop=1000000
c     dtout=1.d-4
c     nstop=1
c----------------------------------------------------------------------|
c  initialize counters
      ns    = 0
      nd=1
      t  = 0.0
      tp = 0.0
      mwflag=0
c----------------------------------------------------------------------|
c   for MPI

      call mpi_init(merrmpi)
      call mpi_comm_size(mpi_comm_world,mpex0,merrmpi)
      call mpi_comm_rank(mpi_comm_world,mpe  ,merrmpi)

      call ix2igx(ix,jx,margin
     &  ,igx,jgx,ig0,ig1,jg0,jg1
     &  ,mpex,ipex,jpex,mpe,ipe,jpe,mpex0,merr)
      if (merr.ne.0) goto 9999
c----------------------------------------------------------------------|
c  file open for "standart output"
      if (mpe.eq.0) then
      mf_out=7
      open(mf_out,file='out.txt',status='replace',iostat=merr)
      if (merr.ne.0) then
        merr=10001
        goto 9999
      endif
      close(mf_out)
      endif
      call mpi_barrier(mpi_comm_world,mermpi)
c----------------------------------------------------------------------|
c  file open
      write(cno,'(i4.4)') mpe
      mf_params=9
      call dacdefparam(mf_params,'params.txt.'//cno)
      mf_t =10
      call dacdef0s(mf_t,'t.dac.'//cno,6)
      mf_pr=21
      call dacdef2s(mf_pr,'pr.dac.'//cno,6,ix,jx)

      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'jx',jx)
      call dacputparami(mf_params,'margin',margin)
      call dacputparami(mf_params,'mpi',1)
      call dacputparami(mf_params,'mpex',mpex)
      call dacputparami(mf_params,'mpe',mpe)
      call dacputparami(mf_params,'ipex',ipex)
      call dacputparami(mf_params,'ipe',ipe)
      call dacputparami(mf_params,'jpex',jpex)
      call dacputparami(mf_params,'jpe',jpe)

c----------------------------------------------------------------------|
c   setup numerical model (grid, initial conditions, etc.)

      call model(ro,pr,gm,rkap0,margin,x,ix,y,jx,mf_params
     &           ,igx,jgx,ipe,jpe)
      mper=(/0,0/)

      call grdrdy(dx,xm,dxm,x,ix)
      call grdrdy(dy,ym,dym,y,jx)
c----------------------------------------------------------------------|
c  read-data

      ndi=1000
      if (mcont.eq.1) then

      mfi_t=60
      call dacopnr0s(mfi_t,'in/t.dac.'//cno,mtype,nx0)
      mfi_pr=71
      call dacopnr2s(mfi_pr,'in/pr.dac.'//cno,mtype,ix0,jx0,nx0)
      do n=1,ndi
        read(mfi_t,end=9900) t
        read(mfi_pr) pr
      enddo
9900  continue

      endif

c----------------------------------------------------------------------|
c     data output
      mf_x=11
      call dacdef1d(mf_x,'x.dac.'//cno,6,ix)
      write(mf_x) x
      mf_y=12
      call dacdef1d(mf_y,'y.dac.'//cno,6,jx)
      write(mf_y) y
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'rkap0',rkap0)
      mf_ro=13
      call dacdef2d(mf_ro,'ro.dac.'//cno,6,ix,jx)
      write(mf_ro) ro

      write(mf_t) t
      write(mf_pr) pr

      if (mpe.eq.0) then
      write(6,913) ns,t,nd
      open(mf_out,file='out.txt',status='old',form='formatted'
     &    ,position='append')
            write(mf_out,913) ns,t,nd
      close(mf_out)
      endif
      nd=nd+1

c======================================================================|
c     time integration
c======================================================================|
1000  continue  
         ns = ns+1
         mwflag=0

c----------------------------------------------------------------------|
c     obtain time spacing

         dt=0.0001d0
         tp = t
         t = t+dt

c----------------------------------------------------------------------|
c     solve conduction equation

         call p_cndsor(ro,pr,mi,err,dt,gm,rkap0
     &                 ,margin,dx,dxm,ix,dy,dym,jx
     &         ,igx,jgx,ipe,jpe,ipex,jpex,mper)

        if (mod(ns,100).eq.0) write(6,925) ns,mi,err

         call exc_1(margin,pr,ix,jx
     &           ,ipe,jpe,ipex,jpex,mper)
         call bnd(margin,pr,ix,jx
     &           ,ipe,jpe,ipex,jpex)

c----------------------------------------------------------------------|
c     data output

         mw=0
         nt1=int(tp/dtout)
         nt2=int(t/dtout)
         if (nt1.lt.nt2) mw=1
         if (mw.ne.0) then
           write(mf_t) t
           write(mf_pr) pr

      if (mpe.eq.0) then
        write(6,913) ns,t,nd
        open(mf_out,file='out.txt',status='old',form='formatted'
     &    ,position='append')
            write(mf_out,913) ns,t,nd
        close(mf_out)
      endif
            nd=nd+1
            mwflag=1
         endif

c----------------------------------------------------------------------|
c     loop test

      if (ns .lt. nstop .and. t .lt. tend) goto 1000
c======================================================================|
c     epilogue
c======================================================================|
9999  continue

c----------------------------------------------------------------------|
c  data output
      if (mwflag.eq.0) then
           write(mf_t) t
           write(mf_pr) pr
      if (mpe.eq.0) then
      write(6,913) ns,t,nd
      open(mf_out,file='out.txt',status='old',form='formatted'
     &    ,position='append')
            write(mf_out,913) ns,t,nd
      close(mf_out)
      endif
      endif

c----------------------------------------------------------------------|
c  ending message
      if(mpe.eq.0)then
      write(6,915) ns,t
      if (merr.eq.0) then
        write(6,*) '  ### normal stop ###'
      else
        write(6,*) '  ### abnormal stop ###'
        write(6,*) '  merr = ',merr
      endif
      open(mf_out,file='out.txt',status='old',form='formatted'
     &    ,position='append')
      write(mf_out,915) ns,t
      if (merr.eq.0) then
        write(mf_out,*) '  ### normal stop ###'
      else
        write(mf_out,*) '  ### abnormal stop ###'
        write(mf_out,*) '  merr = ',merr
      endif
      close(mf_out)
      endif
      call mpi_barrier(mpi_comm_world,mermpi)

      if (merr.ge.10000) then
        call mpi_abort(mpi_comm_world,1,merrmpi)
      else
        call mpi_finalize(merrmpi)
      endif

      stop
      end
