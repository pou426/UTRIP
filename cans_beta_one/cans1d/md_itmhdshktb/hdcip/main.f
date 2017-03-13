c======================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
      parameter (margin=4)
      parameter (ix=1024+2*margin)

      dimension x(ix),xm(ix),dx(ix),dxm(ix)

      dimension ro(ix)
      dimension vx(ix),vy(ix),by(ix)
      dimension bx(ix),bxm(ix)

      dimension vxm(ix),rodx(ix),vxdxm(ix),vydx(ix)

913   format (1x,' write    ','step=',i8,' time=',e10.3,' nd =',i3)
915   format (1x,' stop     ','step=',i8,' time=',e10.3)
c======================================================================|
c     prologue
c======================================================================|
      merr  = 0
      mcont=0
c----------------------------------------------------------------------|
c     set parameters controling finalization and data-output
      tend=0.08d0
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
      mf_ro=20
      call dacdef1s(mf_ro,'ro.dac',6,ix)
      mf_vx=22
      call dacdef1s(mf_vx,'vx.dac',6,ix)
      mf_vy=23
      call dacdef1s(mf_vy,'vy.dac',6,ix)
      mf_by=24
      call dacdef1s(mf_by,'by.dac',6,ix)

      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'margin',margin)
c----------------------------------------------------------------------|
c   setup numerical model (grid, initial conditions, etc.)

      call model(ro,vx,vy,by,bx,bxm,cs2,margin,x,ix,mf_params)

      call grdrdy(dx,xm,dxm,x,ix)
      call ciprdy_mt(vxm,rodx,vxdxm,vydx,ro,vx,vy,dx,dxm,ix)

      call bnd(margin,ro,vx,vy,by,ix)
      call cipbnd(margin,vxm,vy,by,rodx,vxdxm,vydx,ix)
c     floor=1.d-9
c     call chkdav(n_floor,ro,vx,floor,ix)
c----------------------------------------------------------------------|
c  read-data

      ndi=1000
      if (mcont.eq.1) then

      mfi_t=60
      call dacopnr0s(mfi_t,'in/t.dac',mtype,nx0)
      mfi_ro=70
      call dacopnr1s(mfi_ro,'in/ro.dac',mtype,ix0,nx0)
      mfi_vx=72
      call dacopnr1s(mfi_vx,'in/vx.dac',mtype,ix0,nx0)
      mfi_vy=73
      call dacopnr1s(mfi_vy,'in/vy.dac',mtype,ix0,nx0)
      mfi_by=74
      call dacopnr1s(mfi_by,'in/by.dac',mtype,ix0,nx0)
      do n=1,ndi
        read(mfi_t,end=9900) t
        read(mfi_ro) ro
        read(mfi_vx) vx
        read(mfi_vy) vy
        read(mfi_by) by
      enddo
9900  continue

      endif

c----------------------------------------------------------------------|
c     data output
      mf_x=11
      call dacdef1d(mf_x,'x.dac',6,ix)
      write(mf_x) x
      close(mf_x)
      mf_bx=12
      call dacdef1d(mf_bx,'bx.dac',6,ix)
      write(mf_bx) bx
      close(mf_bx)
      call dacputparamd(mf_params,'cs2',cs2)

      write(mf_t) t
      write(mf_ro) ro
      write(mf_vx) vx
      write(mf_vy) vy
      write(mf_by) by
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

         safety=0.9d0
         dtmin=1.d-10
         call cfl_mt(dt,safety,dtmin,merr,cs2,bx,ro,vx,vy,by,dx,ix)
         if (merr.ne.0) goto 9999
         tp = t
         t = t+dt

c----------------------------------------------------------------------|
c     solve hydrodynamic equations

         cvis=0.75d0
         call cip_mt(ro,vx,vy,by,vxm,rodx,vxdxm,vydx
     &    ,bx,bxm,dt,cvis,cs2,dx,dxm,ix)

         call bnd(margin,ro,vx,vy,by,ix)
         call cipbnd(margin,vxm,vy,by,rodx,vxdxm,vydx,ix)

c        floor=1.d-9
c        call chkdav(n_floor,ro,vx,floor,ix)


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
           write(mf_by) by

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
           write(mf_ro) ro
           write(mf_vx) vx
           write(mf_vy) vy
           write(mf_by) by

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
