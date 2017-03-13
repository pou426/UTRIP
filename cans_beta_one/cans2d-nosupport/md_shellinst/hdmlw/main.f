c======================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
      parameter (ix=203,jx=202)

      dimension x(ix),xm(ix),dx(ix),dxm(ix)
      dimension z(jx),zm(jx),dz(jx),dzm(jx)

      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vz(ix,jx)


c======================================================================|
c     prologue
c======================================================================|
      mcont=0
      ndi=1000

c----------------------------------------------------------------------|
c   parameters

      margin=4

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
      mf_vz=24
      call dacdef2s(mf_vz,'vz.dac',6,ix,jx)

      call dacputparamc(mf_params,'comment','cans2d md_shellinst')
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

      tend=5.0
      dtout=0.5
      nstop=1000000
c     dtout=1.d-10
c     nstop=3

c----------------------------------------------------------------------|
c   setup numerical model (grid, initial conditions, etc.)

      call model(ro,pr,vx,vz,gm,margin,x,ix,z,jx
     &    ,mf_params)
      call bnd(margin,ro,pr,vx,vz,ix,jx)
      floor=1.d-9
      call chkdav(n_floor,ro,vx,vz,floor,ix,jx)
      call chkdav(n_floor,pr,vx,vz,floor,ix,jx)

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
      mfi_vz=73
      call dacopnr2s(mfi_vz,'in/vz.dac',mtype,ix0,jx0,nx0)
      do n=1,ndi
        read(mfi_t,end=9900) time
        read(mfi_ro) ro
        read(mfi_pr) pr
        read(mfi_vx) vx
        read(mfi_vz) vz
      enddo
9900  continue

      endif

c----------------------------------------------------------------------|
c     ready

      call grdrdy(dx,xm,dxm,x,ix)
      call grdrdy(dz,zm,dzm,z,jx)


c----------------------------------------------------------------------|
c     data output 
      mf_x=11
      call dacdef1d(mf_x,'x.dac',6,ix)
      write(mf_x) x
      mf_z=12
      call dacdef1d(mf_z,'z.dac',6,jx)
      write(mf_z) z
      call dacputparamd(mf_params,'gm',gm)

      write(mf_t) time
      write(mf_ro) ro
      write(mf_pr) pr
      write(mf_vx) vx
      write(mf_vz) vz
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
         call cfl_h(dt,safety,dtmin,merr,gm,ro,pr,vx,vz,dx,ix,dz,jx)
         if (merr.ne.0) goto 9999
         timep = time
         time = time+dt

c----------------------------------------------------------------------|
c     solve hydrodynamic equations

c                                                      hdmlw - start >>>
         qav=3.d0
         call mlw_h_c(ro,pr,vx,vz,dt,qav,gm,x,xm,dx,dxm,ix,dz,dzm,jx)
c                                                      hdmlw - end   <<<

         call bnd(margin,ro,pr,vx,vz,ix,jx)
         floor=1.d-9
         call chkdav(n_floor,ro,vx,vz,floor,ix,jx)
         call chkdav(n_floor,pr,vx,vz,floor,ix,jx)


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
      write(mf_vz) vz

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
            write(6,913) ns,time,nd
           write(mf_t) time
           write(mf_ro) ro
           write(mf_pr) pr
           write(mf_vx) vx
           write(mf_vz) vz

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
