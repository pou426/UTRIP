c======================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
      parameter (ipex=1,jpex=2,kpex=1,mpex=ipex*jpex*kpex)
      parameter (margin=1)
      parameter (ix=64/ipex+2*margin)
      parameter (jx=64/jpex+2*margin)
      parameter (kx=64/kpex+2*margin)

      dimension x(ix),xm(ix),dx(ix),dxm(ix)
      dimension y(jx),ym(jx),dy(jx),dym(jx)
      dimension z(kx),zm(kx),dz(kx),dzm(kx)

      dimension ro(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension dtarr(2),dtarrg(2)

      include "mpif.h"
      character cno*4
      dimension mper(3)

913   format (1x,' write    ','step=',i8,' time=',e10.3,' nd =',i3)
915   format (1x,' stop     ','step=',i8,' time=',e10.3)
c======================================================================|
c     prologue
c======================================================================|
      merr  = 0
      mcont=0
c----------------------------------------------------------------------|
c     set parameters controling finalization and data-output
      tend=3.0d0
      dtout=0.1d0
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

      call ix2igx(ix,jx,kx,margin
     &  ,igx,jgx,kgx,ig0,ig1,jg0,jg1,kg0,kg1
     &  ,mpex,ipex,jpex,kpex,mpe,ipe,jpe,kpe,mpex0,merr)
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
      call dacdef3s(mf_ro,'ro.dac.'//cno,6,ix,jx,kx)
      mf_vx=22
      call dacdef3s(mf_vx,'vx.dac.'//cno,6,ix,jx,kx)
      mf_vy=23
      call dacdef3s(mf_vy,'vy.dac.'//cno,6,ix,jx,kx)
      mf_vz=24
      call dacdef3s(mf_vz,'vz.dac.'//cno,6,ix,jx,kx)
      mf_bx=25
      call dacdef3s(mf_bx,'bx.dac.'//cno,6,ix,jx,kx)
      mf_by=26
      call dacdef3s(mf_by,'by.dac.'//cno,6,ix,jx,kx)
      mf_bz=27
      call dacdef3s(mf_bz,'bz.dac.'//cno,6,ix,jx,kx)

      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'jx',jx)
      call dacputparami(mf_params,'kx',kx)
      call dacputparami(mf_params,'margin',margin)
      call dacputparami(mf_params,'mpi',1)
      call dacputparami(mf_params,'mpex',mpex)
      call dacputparami(mf_params,'mpe',mpe)
      call dacputparami(mf_params,'ipex',ipex)
      call dacputparami(mf_params,'ipe',ipe)
      call dacputparami(mf_params,'jpex',jpex)
      call dacputparami(mf_params,'jpe',jpe)
      call dacputparami(mf_params,'kpex',kpex)
      call dacputparami(mf_params,'kpe',kpe)

c----------------------------------------------------------------------|
c   setup numerical model (grid, initial conditions, etc.)

      call model(ro,vx,vy,vz,bx,by,bz,cs2,thini,phini
     &           ,margin,x,ix,y,jx,z,kx,mf_params
     &           ,igx,jgx,kgx,ipe,jpe,kpe)
      call pertub(thini,phini,vx,vy,vz,x,ix,y,jx,z,kx)
      mper=(/1,1,1/)
      call grdrdy(dx,xm,dxm,x,ix)
      call grdrdy(dy,ym,dym,y,jx)
      call grdrdy(dz,zm,dzm,z,kx)
c     call exc_7(margin,ro,vx,vy,vz,bx,by,bz,ix,jx,kx
c    &       ,ipe,jpe,kpe,ipex,jpex,kpex,mper)
c     call bnd(margin,ro,vx,vy,vz,bx,by,bz,ix,jx,kx
c    &           ,ipe,jpe,kpe,ipex,jpex,kpex)
c     floor=1.d-9
c     call chkdav(n_floor,ro,vx,vy,vz,floor,ix,jx,kx)

c----------------------------------------------------------------------|
c  read-data

      ndi=1000
      if (mcont.eq.1) then

      mfi_x=61
      call dacopnr1d(mfi_x,'in/x.dac.'//cno,mtype,ix0)
      if (ix.ne.ix0) then
        write(6,*) 'Please change ix :'
        write(6,*) ' ix= ',ix,' ix0= ',ix0
        stop
      endif
      mfi_y=62
      call dacopnr1d(mfi_y,'in/y.dac.'//cno,mtype,jx0)
      if (jx.ne.jx0) then
        write(6,*) 'Please change jx :'
        write(6,*) ' jx= ',jx,' jx0= ',jx0
        stop
      endif
      mfi_z=63
      call dacopnr1d(mfi_z,'in/z.dac.'//cno,mtype,kx0)
      if (kx.ne.kx0) then
        write(6,*) 'Please change kx :'
        write(6,*) ' kx= ',kx,' kx0= ',kx0
        stop
      endif
      mfi_t=60
      call dacopnr0s(mfi_t,'in/t.dac.'//cno,mtype,nx0)
      mfi_ro=70
      call dacopnr3s(mfi_ro,'in/ro.dac.'//cno,mtype,ix0,jx0,kx0,nx0)
      mfi_vx=72
      call dacopnr3s(mfi_vx,'in/vx.dac.'//cno,mtype,ix0,jx0,kx0,nx0)
      mfi_vy=73
      call dacopnr3s(mfi_vy,'in/vy.dac.'//cno,mtype,ix0,jx0,kx0,nx0)
      mfi_vz=74
      call dacopnr3s(mfi_vz,'in/vz.dac.'//cno,mtype,ix0,jx0,kx0,nx0)
      mfi_bx=75
      call dacopnr3s(mfi_bx,'in/bx.dac.'//cno,mtype,ix0,jx0,kx0,nx0)
      mfi_by=76
      call dacopnr3s(mfi_by,'in/by.dac.'//cno,mtype,ix0,jx0,kx0,nx0)
      mfi_bz=77
      call dacopnr3s(mfi_bz,'in/bz.dac.'//cno,mtype,ix0,jx0,kx0,nx0)
      do n=1,ndi
        read(mfi_t,end=9900) t
        read(mfi_ro) ro
        read(mfi_vx) vx
        read(mfi_vy) vy
        read(mfi_vz) vz
        read(mfi_bx) bx
        read(mfi_by) by
        read(mfi_bz) bz
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
      mf_z=13
      call dacdef1d(mf_z,'z.dac.'//cno,6,kx)
      write(mf_z) z
      call dacputparamd(mf_params,'cs2',cs2)

      write(mf_t) t
      write(mf_ro) ro
      write(mf_vx) vx
      write(mf_vy) vy
      write(mf_vz) vz
      write(mf_bx) bx
      write(mf_by) by
      write(mf_bz) bz

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

         safety=0.4d0
         dtmin=1.d-10
         call cfl_mt(dt,safety,dtmin,merr,cs2,ro,vx,vy,vz,bx,by,bz
     &              ,dx,ix,dy,jx,dz,kx)
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

         qav=3.d0
         call mlw_mt(ro,vx,vy,vz,bx,by,bz,dt,qav,cs2
     &            ,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)

      call exc_7(margin,ro,vx,vy,vz,bx,by,bz,ix,jx,kx
     &       ,ipe,jpe,kpe,ipex,jpex,kpex,mper)
      call bnd(margin,ro,vx,vy,vz,bx,by,bz,ix,jx,kx
     &           ,ipe,jpe,kpe,ipex,jpex,kpex)
c        floor=1.d-9
c        call chkdav(n_floor,ro,vx,vy,vz,floor,ix,jx,kx)


c----------------------------------------------------------------------|
c     data output

         mw=0
         nt1=int(tp/dtout)
         nt2=int(t/dtout)
         if (nt1.lt.nt2) mw=1
         if (mw.ne.0) then
           write(mf_t) t
           write(mf_ro) ro
           write(mf_vx) vx
           write(mf_vy) vy
           write(mf_vz) vz
           write(mf_bx) bx
           write(mf_by) by
           write(mf_bz) bz

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
           write(mf_vx) vx
           write(mf_vy) vy
           write(mf_vz) vz
           write(mf_bx) bx
           write(mf_by) by
           write(mf_bz) bz

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
