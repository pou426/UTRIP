c======================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
      parameter (margin=1)
      parameter (ix=128+2*margin,jx=128+2*margin)

      dimension x(ix),xm(ix),dx(ix),dxm(ix)
      dimension z(jx),zm(jx),dz(jx),dzm(jx)

      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx),vz(ix,jx)
      dimension bx(ix,jx),by(ix,jx),bz(ix,jx),ay(ix,jx)

      dimension gx(ix,jx),gxm(ix,jx)
      dimension gz(ix,jx),gzm(ix,jx)

913   format (1x,' write    ','step=',i8,' time=',e10.3,' nd =',i3)
915   format (1x,' stop     ','step=',i8,' time=',e10.3)
c======================================================================|
c     prologue
c======================================================================|
      merr  = 0
      mcont=0
c----------------------------------------------------------------------|
c     set parameters controling finalization and data-output
      tend=6.0d0
      dtout=0.5d0
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
      call dacdef2s(mf_ro,'ro.dac',6,ix,jx)
      mf_pr=21
      call dacdef2s(mf_pr,'pr.dac',6,ix,jx)
      mf_vx=22
      call dacdef2s(mf_vx,'vx.dac',6,ix,jx)
      mf_vy=23
      call dacdef2s(mf_vy,'vy.dac',6,ix,jx)
      mf_vz=24
      call dacdef2s(mf_vz,'vz.dac',6,ix,jx)
      mf_bx=25
      call dacdef2s(mf_bx,'bx.dac',6,ix,jx)
      mf_by=26
      call dacdef2s(mf_by,'by.dac',6,ix,jx)
      mf_bz=27
      call dacdef2s(mf_bz,'bz.dac',6,ix,jx)
      mf_ay=28
      call dacdef2s(mf_ay,'ay.dac',6,ix,jx)

      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'jx',jx)
      call dacputparami(mf_params,'margin',margin)

c----------------------------------------------------------------------|
c   setup numerical model (grid, initial conditions, etc.)

      call model(ro,pr,vx,vy,vz,bx,by,bz,gm
     &     ,gx,gxm,gz,gzm,margin,x,ix,z,jx
     &    ,mf_params)

      call grdrdy(dx,xm,dxm,x,ix)
      call grdrdy(dz,zm,dzm,z,jx)
      call bbtoaa_c(ay,bz,bx,dzm,dxm,x,margin,ix,jx)

      call bnd(margin,ro,pr,vx,vy,vz,bx,by,bz,ay,ix,jx,x,z)
c     floor=1.d-9
c     call chkdav(n_floor,ro,vx,vz,floor,ix,jx)
c     call chkdav(n_floor,pr,vx,vz,floor,ix,jx)
c----------------------------------------------------------------------|
c  read-data

      ndi=1000
      if (mcont.eq.1) then

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
      mfi_vz=74
      call dacopnr2s(mfi_vz,'in/vz.dac',mtype,ix0,jx0,nx0)
      mfi_bx=75
      call dacopnr2s(mfi_bx,'in/bx.dac',mtype,ix0,jx0,nx0)
      mfi_by=76
      call dacopnr2s(mfi_by,'in/by.dac',mtype,ix0,jx0,nx0)
      mfi_bz=77
      call dacopnr2s(mfi_bz,'in/bz.dac',mtype,ix0,jx0,nx0)
      mfi_ay=78
      call dacopnr2s(mfi_ay,'in/ay.dac',mtype,ix0,jx0,nx0)
      do n=1,ndi
        read(mfi_t,end=9900) t
        read(mfi_ro) ro
        read(mfi_pr) pr
        read(mfi_vx) vx
        read(mfi_vy) vy
        read(mfi_vz) vz
        read(mfi_bx) bx
        read(mfi_by) by
        read(mfi_bz) bz
        read(mfi_ay) ay
      enddo
9900  continue

      endif

c----------------------------------------------------------------------|
c     data output
      mf_x=11
      call dacdef1d(mf_x,'x.dac',6,ix)
      write(mf_x) x
      mf_z=12
      call dacdef1d(mf_z,'z.dac',6,jx)
      write(mf_z) z
      call dacputparamd(mf_params,'gm',gm)

      write(mf_t) t
      write(mf_ro) ro
      write(mf_pr) pr
      write(mf_vx) vx
      write(mf_vy) vy
      write(mf_vz) vz
      write(mf_bx) bx
      write(mf_by) by
      write(mf_bz) bz
      write(mf_ay) ay
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

         safety=0.4d0
         dtmin=1.d-10
         call cfl_m3(dt,safety,dtmin,merr,gm,ro,pr,vx,vy,vz,bx,by,bz
     &      ,dx,ix,dz,jx)
         if (merr.ne.0) goto 9999
         tp = t
         t = t+dt

c----------------------------------------------------------------------|
c     solve hydrodynamic equations

         qav=3.d0
         call mlw_m3_cg(ro,pr,vx,vy,vz,bx,by,bz,ay,dt,qav,gm
     &                ,gx,gxm,gz,gzm,x,xm,dx,dxm,ix,dz,dzm,jx)

         call bnd(margin,ro,pr,vx,vy,vz,bx,by,bz,ay,ix,jx,x,z)
c        floor=1.d-9
c        call chkdav(n_floor,ro,vx,vz,floor,ix,jx)
c        call chkdav(n_floor,pr,vx,vz,floor,ix,jx)


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
           write(mf_vy) vy
           write(mf_vz) vz
           write(mf_bx) bx
           write(mf_by) by
           write(mf_bz) bz
           write(mf_ay) ay

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
           write(mf_pr) pr
           write(mf_vx) vx
           write(mf_vy) vy
           write(mf_vz) vz
           write(mf_bx) bx
           write(mf_by) by
           write(mf_bz) bz
           write(mf_ay) ay

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
