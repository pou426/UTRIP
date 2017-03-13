c======================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
      parameter (margin=1)
      parameter (ix=1024+2*margin)

      dimension x(ix),xm(ix),dx(ix),dxm(ix)
      dimension ey(ix),ez(ix),by(ix),bz(ix)

913   format (1x,' write    ','step=',i8,' time=',e10.3,' nd =',i3)
915   format (1x,' stop     ','step=',i8,' time=',e10.3)
c======================================================================|
c     prologue
c======================================================================|
      merr  = 0
      mcont=0
c----------------------------------------------------------------------|
c     set parameters controling finalization and data-output
      tend=1.0d0
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
c  file open for "standart output"
      mf_out=7
      open(mf_out,file='out.txt',status='replace',iostat=merr)
      if (merr.ne.0) then
        merr=10001
        goto 9999
      endif
      close(mf_out)
c----------------------------------------------------------------------|
c  file open
      mf_params=9
      call dacdefparam(mf_params,'params.txt')
      mf_t =10
      call dacdef0s(mf_t,'t.dac',6)
      mf_ey=20
      call dacdef1s(mf_ey,'ey.dac',6,ix)
      mf_ez=21
      call dacdef1s(mf_ez,'ez.dac',6,ix)
      mf_by=22
      call dacdef1s(mf_by,'by.dac',6,ix)
      mf_bz=23
      call dacdef1s(mf_bz,'bz.dac',6,ix)

      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'margin',margin)
c----------------------------------------------------------------------|
c   setup numerical model (grid, initial conditions, etc.)

      call model(ey,ez,by,bz,margin,x,ix,mf_params)

      call grdrdy(dx,xm,dxm,x,ix)

      call bnd(margin,ey,ez,by,bz,ix)
c----------------------------------------------------------------------|
c  read-data

      ndi=1000
      if (mcont.eq.1) then

      mfi_t=60
      call dacopnr0s(mfi_t,'in/t.dac',mtype,nx0)
      mfi_ey=70
      call dacopnr1s(mfi_ey,'in/ey.dac',mtype,ix0,nx0)
      mfi_ez=71
      call dacopnr1s(mfi_ez,'in/ez.dac',mtype,ix0,nx0)
      mfi_by=72
      call dacopnr1s(mfi_by,'in/by.dac',mtype,ix0,nx0)
      mfi_bz=73
      call dacopnr1s(mfi_bz,'in/bz.dac',mtype,ix0,nx0)
      do n=1,ndi
        read(mfi_t,end=9900) t
        read(mfi_ey) ey
        read(mfi_ez) ez
        read(mfi_by) by
        read(mfi_bz) bz
      enddo
9900  continue

      endif

c----------------------------------------------------------------------|
c     data output
      mf_x=11
      call dacdef1d(mf_x,'x.dac',6,ix)
      write(mf_x) x
      close(mf_x)

      write(mf_t) t
      write(mf_ey) ey
      write(mf_ez) ez
      write(mf_by) by
      write(mf_bz) bz
      write(6,913) ns,t,nd
      open(mf_out,file='out.txt',status='old',form='formatted'
     &    ,position='append')
            write(mf_out,913) ns,t,nd
      close(mf_out)
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
         call cfl_em(dt,safety,dtmin,merr,dx,ix)
         if (merr.ne.0) goto 9999
         tp = t
         t = t+dt

c----------------------------------------------------------------------|
c     solve hydrodynamic equations

         call mlw_em(ey,ez,by,bz,dt,dx,dxm,ix)

         call bnd(margin,ey,ez,by,bz,ix)

c----------------------------------------------------------------------|
c     data output

         mw=0
         nt1=int(tp/dtout)
         nt2=int(t/dtout)
         if (nt1.lt.nt2) mw=1
         if (mw.ne.0) then
           write(mf_t) t
            write(mf_ey) ey
            write(mf_ez) ez
            write(mf_by) by
            write(mf_bz) bz

           write(6,913) ns,t,nd
           open(mf_out,file='out.txt',status='old',form='formatted'
     &         ,position='append')
                 write(mf_out,913) ns,t,nd
           close(mf_out)
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
            write(mf_ey) ey
            write(mf_ez) ez
            write(mf_by) by
            write(mf_bz) bz

           write(6,913) ns,t,nd
           open(mf_out,file='out.txt',status='old',form='formatted'
     &         ,position='append')
                 write(mf_out,913) ns,t,nd
           close(mf_out)
      endif

c----------------------------------------------------------------------|
c  ending message
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

      stop
      end
