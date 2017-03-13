c=====================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
      parameter (margin=1)
      parameter (ix=1024+2*margin)

      dimension x(ix),xm(ix),dx(ix),dxm(ix)

      dimension ro(ix),pr(ix)

      dimension sc(ix),scm(ix)

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
      tend=10.0d0
      dtout=1.0d0
      nstop=1000000
c     dtout=1.d-4
c     nstop=1
c----------------------------------------------------------------------|
c  initialize counters
      ns    = 0
      t  = 0.0
      tp = 0.0
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
      mf_pr=20
      call dacdef1s(mf_pr,'pr.dac',6,ix)

      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'margin',margin)

c----------------------------------------------------------------------|
c   setup numerical model (grid, initial conditions, etc.)

      call model(ro,pr,gm,rkap0,margin,sc,scm,x,ix,mf_params)
      
      call grdrdy(dx,xm,dxm,x,ix)
c     call bnd(margin,pr,ix)

c----------------------------------------------------------------------|
c  read-data

      ndi=1000
      if (mcont.eq.1) then

      mfi_t=60
      call dacopnr0s(mfi_t,'in/t.dac',mtype,nx0)
      mfi_pr=70
      call dacopnr1s(mfi_pr,'in/pr.dac',mtype,ix0,nx0)
      do n=1,ndi
        read(mfi_t,end=9900) t
        read(mfi_pr) pr
      enddo
9900  continue

      endif

c----------------------------------------------------------------------|
c     data output
      mf_x=11
      call dacdef1d(mf_x,'x.dac',6,ix)
      write(mf_x) x
      close(mf_x)
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'rkap0',rkap0)
      mf_ro=12
      call dacdef1d(mf_ro,'ro.dac',6,ix)
      write(mf_ro) ro
      close(mf_ro)
      mf_sc=13
      call dacdef1d(mf_sc,'sc.dac',6,ix)
      write(mf_sc) sc
      close(mf_sc)
      mf_scm=14
      call dacdef1d(mf_scm,'scm.dac',6,ix)
      write(mf_scm) scm
      close(mf_scm)

      write(mf_t) t
      write(mf_pr) pr

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

         dt=0.01d0
         tp = t
         t = t+dt

c----------------------------------------------------------------------|
c     solve conduction equation

         call cndsor_c(ro,pr,mi,err,dt,gm,rkap0,margin
     &          ,sc,scm,dx,dxm,ix)

        if (mod(ns,1000).eq.0) write(6,925) ns,mi,err

         call bnd(margin,pr,ix)

c----------------------------------------------------------------------|
c     data output

         mw=0
         nt1=int(tp/dtout)
         nt2=int(t/dtout)
         if (nt1.lt.nt2) mw=1
         if (mw.ne.0) then
           write(mf_t) t
           write(mf_pr) pr

            write(6,913) ns,t,nd
           open(mf_out,file='out.txt',status='old',form='formatted'
     &    ,position='append')
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
           write(mf_pr) pr
            write(6,913) ns,t,nd
      open(mf_out,file='out.txt',status='old',form='formatted'
     &    ,position='append')
            write(mf_out,913) ns,t,nd
      close(mf_out)

      endif

c----------------------------------------------------------------------|
c  file close

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