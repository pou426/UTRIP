c======================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
c     parameter (ix=128,jx=3)
c     parameter (ix=3,jx=128)
      parameter (ix=128,jx=128)

      dimension x(ix),xm(ix),dx(ix),dxm(ix)
      dimension y(jx),ym(jx),dy(jx),dym(jx)

      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx)
      dimension bx(ix,jx),by(ix,jx)
      dimension az(ix,jx)


c======================================================================|
c     prologue
c======================================================================|
      mcont=0
      ndi=1000

c----------------------------------------------------------------------|
c   parameters

      margin=1

c----------------------------------------------------------------------|
c  file open
      mf_params=9
      call dacdefparam(mf_params,'params.txt')
      mf_t =10
      call dacdef0s(mf_t,'t.dac',6)
      mf_ro=20
      call dacdef2s(mf_ro,'ro.dac',6,ix,jx)
      mf_pr=21
      call dacdef2s(mf_pr,'pr.dac',6,ix,jx)
      mf_vx=22
      call dacdef2s(mf_vx,'vx.dac',6,ix,jx)
      mf_vy=23
      call dacdef2s(mf_vy,'vy.dac',6,ix,jx)
      mf_bx=25
      call dacdef2s(mf_bx,'bx.dac',6,ix,jx)
      mf_by=26
      call dacdef2s(mf_by,'by.dac',6,ix,jx)
      mf_az=28
      call dacdef2s(mf_az,'az.dac',6,ix,jx)

      call dacputparamc(mf_params,'comment','cans2d md_relmhdshk')
      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'jx',jx)
      call dacputparami(mf_params,'margin',margin)

c----------------------------------------------------------------------|
c  initialize counters

      nd=1
      time  = 0.0
      timep = 0.0
      ns    = 0
      merr  = 0

c----------------------------------------------------------------------|
c   time control parameters

      dtout=1.d0
      tend=10.d0
      nstop=100000
c     dtout=1.e-6
c     nstop=1
c----------------------------------------------------------------------|
c   setup numerical model (grid, initial conditions, etc.)

      call model(ro,pr,vx,vy,bx,by,gm,margin,x,ix,y,jx
     &    ,tend,dtout,mf_params)

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
      mfi_y=62
      call dacopnr1d(mfi_y,'in/y.dac',mtype,jx0)
      if (jx.ne.jx0) then
        write(6,*) 'Please change jx :'
        write(6,*) ' jx= ',jx,' jx0= ',jx0
        stop
      endif
      mfi_t=60
      call dacopnr0s(mfi_t,'in/t.dac',mtype,nx0)
      mfi_ro=70
      call dacopnr2s(mfi_ro,'in/ro.dac',mtype,ix0,jx0,nx0)
      mfi_pr=71
      call dacopnr2s(mfi_pr,'in/pr.dac',mtype,ix0,jx0,nx0)
      mfi_vx=72
      call dacopnr2s(mfi_vx,'in/vx.dac',mtype,ix0,jx0,nx0)
      mfi_vy=73
      call dacopnr2s(mfi_vy,'in/vy.dac',mtype,ix0,jx0,nx0)
      mfi_bx=74
      call dacopnr2s(mfi_bx,'in/bx.dac',mtype,ix0,jx0,nx0)
      mfi_by=75
      call dacopnr2s(mfi_by,'in/by.dac',mtype,ix0,jx0,nx0)
      mfi_az=76
      call dacopnr2s(mfi_az,'in/az.dac',mtype,ix0,jx0,nx0)
      do n=1,ndi
        read(mfi_t,end=9900) time
        read(mfi_ro) ro
        read(mfi_pr) pr
        read(mfi_vx) vx
        read(mfi_vy) vy
        read(mfi_bx) bx
        read(mfi_by) by
        read(mfi_az) az
      enddo
9900  continue

      endif

c----------------------------------------------------------------------|
c     ready

      call grdrdy(dx,xm,dxm,x,ix)
      call grdrdy(dy,ym,dym,y,jx)
      call bbtoaa(az,bx,by,dxm,dym,ix,jx)
      call bnd(margin,ro,pr,vx,vy,bx,by,az,ix,jx)
      floor=1.d-9
      call chkdav(n_floor,ro,vx,vy,floor,ix,jx)
      call chkdav(n_floor,pr,vx,vy,floor,ix,jx)


c----------------------------------------------------------------------|
c     data output
      mf_x=11
      call dacdef1d(mf_x,'x.dac',6,ix)
      write(mf_x) x
      mf_y=12
      call dacdef1d(mf_y,'y.dac',6,jx)
      write(mf_y) y
      call dacputparamd(mf_params,'gm',gm)

      write(mf_t) time
      write(mf_ro) ro
      write(mf_pr) pr
      write(mf_vx) vx
      write(mf_vy) vy
      write(mf_bx) bx
      write(mf_by) by
      write(mf_az) az
      write(6,913) ns,time,nd
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
         call cfl_m(dt,safety,dtmin,merr,gm,ro,pr,vx,vy,bx,by
     &      ,dx,ix,dy,jx)

         if (merr.ne.0) goto 9999
         timep = time
         time = time+dt

c----------------------------------------------------------------------|
c     solve hydrodynamic equations

c                                                      hdmlw - start >>>
         qav=3.d0
         call mlw_sm(ro,pr,vx,vy,bx,by,az,dt,qav,gm,dx,dxm,ix,dy,dym,jx)
c                                                      hdmlw - end   <<<

         call bnd(margin,ro,pr,vx,vy,bx,by,az,ix,jx)
         floor=1.d-9
         call chkdav(n_floor,ro,vx,vy,floor,ix,jx)
         call chkdav(n_floor,pr,vx,vy,floor,ix,jx)


c----------------------------------------------------------------------|
c     data output

         mw=0
         nt1=int(timep/dtout)
         nt2=int(time/dtout)
         if (nt1.lt.nt2) mw=1
         if (mw.ne.0) then
      write(mf_t) time
      write(mf_ro) ro
      write(mf_pr) pr
      write(mf_vx) vx
      write(mf_vy) vy
      write(mf_bx) bx
      write(mf_by) by
      write(mf_az) az

            write(6,913) ns,time,nd
            nd=nd+1
            mwflag=1
         endif

c----------------------------------------------------------------------|
c     loop test

      if (ns .lt. nstop .and. time .lt. tend) goto 1000
         
c======================================================================|
c     epilogue
c======================================================================|
9999  continue

c----------------------------------------------------------------------|
c  data output
      if (mwflag.eq.0) then
           write(mf_t) time
           write(mf_ro) ro
           write(mf_pr) pr
           write(mf_vx) vx
           write(mf_vy) vy
           write(mf_bx) bx
           write(mf_by) by
           write(mf_az) az

            write(6,913) ns,time,nd
      endif

c----------------------------------------------------------------------|
c  file close

c----------------------------------------------------------------------|
c  ending message

      write(6,915) ns,time
      if (merr.eq.0) then
        write(6,*) '  ### normal stop ###'
      else
        write(6,*) '  ### abnormal stop ###'
      endif

      stop

913   format (1x,' write    ','step=',i8,' time=',e10.3,' nd =',i3)
915   format (1x,' stop     ','step=',i8,' time=',e10.3)

      end
