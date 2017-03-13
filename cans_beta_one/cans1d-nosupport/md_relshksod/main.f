c======================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
      parameter (ix=508)

      dimension x(ix),xm(ix),dx(ix),dxm(ix)
      dimension ro(ix),pr(ix),vx(ix)


c                                                      hdcip - start >>>
c     dimension de(ix),ei(ix),rxm(ix),dedx(ix),eidx(ix),rxdxm(ix)
c                                                      hdcip - end   <<<

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

      dtout=0.05d0
      tend=0.14154d0
      dtout=1.d0
      tend=45.d0
      dtout=500.d0
      tend=5000.d0
      nstop=50000
c     dtout=1.d-6
c     nstop=10

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

      call dacputparamc(mf_params,'comment','cans1d md_relshksod')
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

      call model(ro,pr,vx,gm,margin,x,ix
     &     ,mf_params)
      call bnd(margin,ro,pr,vx,ix)

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
      do n=1,ndi
        read(mfi_t,end=9900) t
        read(mfi_ro) ro
        read(mfi_pr) pr
        read(mfi_vx) vx
      enddo
9900  continue
      close(mfi_t)
      close(mfi_ro)
      close(mfi_pr)
      close(mfi_vx)

      endif

c----------------------------------------------------------------------|
c     ready

      call grdrdy(dx,xm,dxm,x,ix)
c                                                      hdcip - start >>>
c     call cippre_sh(ro,pr,vx,de,ei,rxm,gm,ix)
c     call ciprdy_sh(de,ei,rxm,dedx,eidx,rxdxm,dx,dxm,ix)
c     call cipbnd(margin,de,ei,rxm,dedx,eidx,rxdxm,ix)
c                                                      hdcip - end   <<<

c----------------------------------------------------------------------|
c     data output 
      mf_x=11
      call dacdef1d(mf_x,'x.dac',6,ix)
      write(mf_x) x
      close(mf_x)
      call dacputparamd(mf_params,'gm',gm)

      write(mf_t) t
      write(mf_ro) ro
      write(mf_pr) pr
      write(mf_vx) vx

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
         call cfl_h(dt,safety,dtmin,merr,gm,ro,pr,vx,dx,ix)
         if (merr.ne.0) goto 9999
         tp = t
         t = t+dt

c----------------------------------------------------------------------|
c     solve hydrodynamic equations

c                                                      hdcip - start >>>
c        call cippre_sh(ro,pr,vx,de,ei,rxm,gm,ix)
c        cvis=0.75d0
c        call cip_sh(de,ei,rxm,dedx,eidx,rxdxm,dt,cvis,gm,dx,dxm,ix)
c        call cippost_sh(de,ei,rxm,ro,pr,vx,gm,ix)
c        call cipbnd(margin,de,ei,rxm,dedx,eidx,rxdxm,ix)
c                                                      hdcip - end   <<<

         qav=3.d0
         call mlw_sh(ro,pr,vx,dt,qav,gm,dx,dxm,ix)

         call bnd(margin,ro,pr,vx,ix)
         floor=1.d-9
         call chkdav(n_floor,ro,vx,floor,ix)
         call chkdav(n_floor,pr,vx,floor,ix)

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
      endif

c----------------------------------------------------------------------|
c  file close
      close(mf_t)
      close(mf_ro)
      close(mf_pr)
      close(mf_vx)

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
