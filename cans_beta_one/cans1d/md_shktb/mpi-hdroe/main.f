c======================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
      parameter (ipex=2,mpex=ipex)
      parameter (margin=3)
      parameter (ix=1024/ipex+2*margin)

      dimension x(ix),xm(ix),dx(ix),dxm(ix)

      dimension ro(ix),pr(ix)
      dimension vx(ix)
      dimension dtarr(2),dtarrg(2)

      include "mpif.h"
      character cno*4

913   format (1x,' write    ','step=',i8,' time=',e10.3,' nd =',i3)
915   format (1x,' stop     ','step=',i8,' time=',e10.3)
c======================================================================|
c     prologue
c======================================================================|
      merr  = 0
      mcont=0
c----------------------------------------------------------------------|
c     set parameters controling finalization and data-output
      tend=0.14154d0
      dtout=0.02d0
      nstop=100000
c     dtout=1.d-10
c     nstop=1
c----------------------------------------------------------------------|
c  initialize counters
      ns    = 0
      t  = 0.0d0
      tp = 0.0d0
      nd=1
      mwflag=0
c----------------------------------------------------------------------|
c   for MPI

      call mpi_init(merrmpi)
      call mpi_comm_size(mpi_comm_world,mpex0,merrmpi)
      call mpi_comm_rank(mpi_comm_world,mpe  ,merrmpi)

      call ix2igx(ix,margin
     &  ,igx,ig0,ig1
     &  ,mpex,ipex,mpe,ipe,mpex0,merr)
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
      mf_ro=20
      call dacdef1s(mf_ro,'ro.dac.'//cno,6,ix)
      mf_pr=21
      call dacdef1s(mf_pr,'pr.dac.'//cno,6,ix)
      mf_vx=22
      call dacdef1s(mf_vx,'vx.dac.'//cno,6,ix)

      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'margin',margin)
      call dacputparami(mf_params,'mpi',1)
      call dacputparami(mf_params,'mpex',mpex)
      call dacputparami(mf_params,'mpe',mpe)
      call dacputparami(mf_params,'ipex',ipex)
      call dacputparami(mf_params,'ipe',ipe)

c----------------------------------------------------------------------|
c   setup numerical model (grid, initial conditions, etc.)

      call model(ro,pr,vx,gm,margin,x,dxm,ix,mf_params,igx,ipe)
      mper=0

      call grdrdy(dx,xm,dxm,x,ix)
c     call exc_3(margin,ro,pr,vx,ix,ipe,ipex,mper)
c     call bnd(margin,ro,pr,vx,ix,ipe,ipex)
c     floor=1.d-9
c     call chkdav(n_floor,ro,vx,floor,ix)
c     call chkdav(n_floor,pr,vx,floor,ix)

c----------------------------------------------------------------------|
c  read-data

      ndi=1000
      if (mcont.eq.1) then

      mfi_t=60
      call dacopnr0s(mfi_t,'in/t.dac.'//cno,mtype,nx0)
      mfi_ro=70
      call dacopnr1s(mfi_ro,'in/ro.dac.'//cno,mtype,ix0,nx0)
      mfi_pr=71
      call dacopnr1s(mfi_pr,'in/pr.dac.'//cno,mtype,ix0,nx0)
      mfi_vx=72
      call dacopnr1s(mfi_vx,'in/vx.dac.'//cno,mtype,ix0,nx0)
      do n=1,ndi
        read(mfi_t,end=9900) t
        read(mfi_ro) ro
        read(mfi_pr) pr
        read(mfi_vx) vx
      enddo
9900  continue

      endif

c----------------------------------------------------------------------|
c     data output
      mf_x=11
      call dacdef1d(mf_x,'x.dac.'//cno,6,ix)
      write(mf_x) x
      close(mf_x)
      call dacputparamd(mf_params,'gm',gm)

      write(mf_t) t
      write(mf_ro) ro
      write(mf_pr) pr
      write(mf_vx) vx

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

         safety=1.0d0
         dtmin=1.d-10
         call cfl_h(dt,safety,dtmin,merr,gm,ro,pr,vx,dx,ix)
         dtarr=(/dt,dble(-merr)/)
         dtarrg=(/0.d0,0.d0/)
         call mpi_allreduce(dtarr,dtarrg,2,mpi_double_precision,mpi_min
     &                      ,mpi_comm_world,merrmpi)
         merr = int(-dtarrg(2)*1.01)
         if (merr.ne.0) goto 9999
         dt = dtarrg(1)
         tp = t
         t  = t+dt

c----------------------------------------------------------------------|
c     solve hydrodynamic equations

         call roe_h(ro,pr,vx,dt,gm,dx,ix)

         call exc_3(margin,ro,pr,vx,ix,ipe,ipex,mper)
         call bnd(margin,ro,pr,vx,ix,ipe,ipex)
c        floor=1.d-9
c        call chkdav(n_floor,ro,vx,floor,ix)
c        call chkdav(n_floor,pr,vx,floor,ix)

c----------------------------------------------------------------------|
c     data output

         mw=0
         nt1=int(tp/dtout)
         nt2=int(t/dtout)
         if (nt1.lt.nt2) mw=1
         if (mw.ne.0) then
           write(mf_t) t
           write(mf_ro) ro
           write(mf_pr) pr
           write(mf_vx) vx

           if (mpe.eq.0) then
             write(6,913) ns,t,nd
             open(mf_out,file='out.txt',status='old',form='formatted'
     &        ,position='append')
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
           write(mf_ro) ro
           write(mf_pr) pr
           write(mf_vx) vx

           if (mpe.eq.0) then
             write(6,913) ns,t,nd
             open(mf_out,file='out.txt',status='old',form='formatted'
     &        ,position='append')
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
