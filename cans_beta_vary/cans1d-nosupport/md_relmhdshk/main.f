c======================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
      parameter (ix=256)

      dimension x(ix),xm(ix),dx(ix),dxm(ix)
      dimension ro(ix),pr(ix),vx(ix),vy(ix),by(ix)
      dimension bx(ix),bxm(ix)

c======================================================================|
c     prologue
c======================================================================|
      mcont=0
      ndi=1000

c----------------------------------------------------------------------|
c   parameters

      margin=4

c----------------------------------------------------------------------|
c   time control parameters

c     dtout=0.1d0
c     tend=2.d0
      nstop=50000
c     dtout=1.d-6
c     nstop=1

c----------------------------------------------------------------------|
c  file open
      mf_params=9
      call dacdefparam(mf_params,'params.txt')
      mf_t =10
      call dacdef0s(mf_t,'t.dac',6)
      mf_ro=20
      call dacdef1s(mf_ro,'ro.dac',6,ix)
      mf_pr=21
      call dacdef1s(mf_pr,'pr.dac',6,ix)
      mf_vx=22
      call dacdef1s(mf_vx,'vx.dac',6,ix)
      mf_vy=23
      call dacdef1s(mf_vy,'vy.dac',6,ix)
      mf_by=24
      call dacdef1s(mf_by,'by.dac',6,ix)

      call dacputparamc(mf_params,'comment','cans1d md_relmhdshk')
      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'margin',margin)
c----------------------------------------------------------------------|
c  initialize counters

      nd=1
      t  = 0.0
      tp = 0.0
      ns    = 0
      merr  = 0

c----------------------------------------------------------------------|
c   setup numerical model (grid, initial conditions, etc.)

      call model(ro,pr,vx,vy,by,bx,bxm,gm,margin,x,ix
     &     ,tend,dtout,mf_params)
      call bnd(margin,ro,pr,vx,vy,by,ix)

c----------------------------------------------------------------------|
c  read-data

      if (mcont.eq.1) then

      mfi_x=61
      call dacopnr1d(mfi_x,'in/x.dac',mtype,ix0)
      if (ix.ne.ix0) then
        write(6,*) 'Please change ix :'
        write(6,*) ' ix= ',ix,' ix0= ',ix0
        stop
      endif
      mfi_t=60
      call dacopnr0s(mfi_t,'in/t.dac',mtype,nx0)
      mfi_ro=70
      call dacopnr1s(mfi_ro,'in/ro.dac',mtype,ix0,nx0)
      mfi_pr=71
      call dacopnr1s(mfi_pr,'in/pr.dac',mtype,ix0,nx0)
      mfi_vx=72
      call dacopnr1s(mfi_vx,'in/vx.dac',mtype,ix0,nx0)
      mfi_vy=73
      call dacopnr1s(mfi_vy,'in/vy.dac',mtype,ix0,nx0)
      mfi_by=74
      call dacopnr1s(mfi_by,'in/by.dac',mtype,ix0,nx0)
      do n=1,ndi
        read(mfi_t,end=9900) t
        read(mfi_ro) ro
        read(mfi_pr) pr
        read(mfi_vx) vx
        read(mfi_vy) vy
        read(mfi_by) by
      enddo
9900  continue
      close(mfi_t)
      close(mfi_ro)
      close(mfi_pr)
      close(mfi_vx)
      close(mfi_vy)
      close(mfi_by)

      endif

c----------------------------------------------------------------------|
c     ready

      call grdrdy(dx,xm,dxm,x,ix)

c----------------------------------------------------------------------|
c     data output 
      mf_x=11
      call dacdef1d(mf_x,'x.dac',6,ix)
      write(mf_x) x
      close(mf_x)
      call dacputparamd(mf_params,'gm',gm)
      mf_bx=12
      call dacdef1d(mf_bx,'bx.dac',6,ix)
      write(mf_bx) bx
      close(mf_bx)

      write(mf_t) t
      write(mf_ro) ro
      write(mf_pr) pr
      write(mf_vx) vx
      write(mf_vy) vy
      write(mf_by) by

      write(6,913) ns,t,nd
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
         call cfl_m(dt,safety,dtmin,merr,gm,bx,ro,pr,vx,vy,by,dx,ix)
         if (merr.ne.0) goto 9999
         tp = t
         t = t+dt

c----------------------------------------------------------------------|
c     solve hydrodynamic equations

c                                                      hdmlw - start >>>
         qav=3.d0
         call mlw_sm(ro,pr,vx,vy,by,bx,bxm,dt,qav,gm,dx,dxm,ix,ns)
c                                                      hdmlw - end   <<<
         call bnd(margin,ro,pr,vx,vy,by,ix)
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
            write(6,913) ns,t,nd
            write(mf_t) t
            write(mf_ro) ro
            write(mf_pr) pr
            write(mf_vx) vx
            write(mf_vy) vy
            write(mf_by) by

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
          write(6,913) ns,t,nd
            write(mf_t) t
            write(mf_ro) ro
            write(mf_pr) pr
            write(mf_vx) vx
            write(mf_vy) vy
            write(mf_by) by
      endif

c----------------------------------------------------------------------|
c  file close
      close(mf_t)
      close(mf_ro)
      close(mf_pr)
      close(mf_vx)
      close(mf_vy)
      close(mf_by)

c----------------------------------------------------------------------|
c  ending message

      write(6,915) ns,t
      if (merr.eq.0) then
        write(6,*) '  ### normal stop ###'
      else
        write(6,*) '  ### abnormal stop ###'
      endif 

      stop

913   format (1x,' write    ','step=',i8,' t=',e10.3,' nd =',i3)
915   format (1x,' stop     ','step=',i8,' t=',e10.3)

      end
