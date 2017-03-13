c======================================================================|
c     array definitions
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension mfdim(3)
      data mfdim/1,1,1/
      parameter (margin=1)
      parameter (ix=64+2*margin,jx=4+2*margin,kx=64+2*margin)

      dimension r(ix),rm(ix),dr(ix),drm(ix)
      dimension ph(jx),phm(jx),dph(jx),dphm(jx)
      dimension z(kx),zm(kx),dz(kx),dzm(kx)

      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vr(ix,jx,kx),vph(ix,jx,kx),vz(ix,jx,kx)
      dimension br(ix,jx,kx),bph(ix,jx,kx),bz(ix,jx,kx)
      double precision mu

913   format (1x,' write    ','step=',i8,' time=',e10.3,' nd =',i3)
915   format (1x,' stop     ','step=',i8,' time=',e10.3)
c======================================================================|
c     prologue
c======================================================================|
      merr  = 0
      mcont=0
c----------------------------------------------------------------------|
c     set parameters controling finalization and data-output
      tend=5.0d0
      dtout=0.5d0
      nstop=100000
c     dtout=1.d-10
c     nstop=1
c----------------------------------------------------------------------|
c  initialize counters
      ns    = 0
      t  = 0.0d0
      tp = 0.0d0
      nd=0
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
      call dacdef3s(mf_ro,'ro.dac',6,ix,jx,kx)
      mf_pr=21
      call dacdef3s(mf_pr,'pr.dac',6,ix,jx,kx)
      mf_vr=22
      call dacdef3s(mf_vr,'vr.dac',6,ix,jx,kx)
      mf_vph=23
      call dacdef3s(mf_vph,'vph.dac',6,ix,jx,kx)
      mf_vz=24
      call dacdef3s(mf_vz,'vz.dac',6,ix,jx,kx)
      mf_br=25
      call dacdef3s(mf_br,'br.dac',6,ix,jx,kx)
      mf_bph=26
      call dacdef3s(mf_bph,'bph.dac',6,ix,jx,kx)
      mf_bz=27
      call dacdef3s(mf_bz,'bz.dac',6,ix,jx,kx)

      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'jx',jx)
      call dacputparami(mf_params,'kx',kx)
      call dacputparami(mf_params,'margin',margin)


c----------------------------------------------------------------------|
c   setup numerical model (grid, initial conditions, etc.)

      call model(ro,pr,vr,vph,vz,br,bph,bz,gm,mu
     &   ,r,dr,rm,drm,ph,dph,phm,dphm,z,dz,zm,dzm
     &   ,margin,ix,jx,kx,mf_params)
c     call grdrdy(dr,rm,drm,r,ix)
c     call grdrdy(dph,phm,dphm,ph,jx)
c     call grdrdy(dz,zm,dzm,z,kx)
c     call bnd(margin,ro,pr,vr,vph,vz,br,bph,bz,ix,jx,kx,mfdim)
c     floor=1.d-9
c     call chkdav(n_floor,ro,vr,vph,vz,floor,ix,jx,kx)
c     call chkdav(n_floor,pr,vr,vph,vz,floor,ix,jx,kx)


c----------------------------------------------------------------------|
c  read-data

      ndi=1000
      if (mcont.eq.1) then

      mfi_t=60
      call dacopnr0s(mfi_t,'in/t.dac',mtype,nx0)
      mfi_ro=70
      call dacopnr3s(mfi_ro,'in/ro.dac',mtype,ix0,jx0,kx0,nx0)
      mfi_pr=71
      call dacopnr3s(mfi_pr,'in/pr.dac',mtype,ix0,jx0,kx0,nx0)
      mfi_vr=72
      call dacopnr3s(mfi_vr,'in/vr.dac',mtype,ix0,jx0,kx0,nx0)
      mfi_vph=73
      call dacopnr3s(mfi_vph,'in/vph.dac',mtype,ix0,jx0,kx0,nx0)
      mfi_vz=74
      call dacopnr3s(mfi_vz,'in/vz.dac',mtype,ix0,jx0,kx0,nx0)
      mfi_br=75
      call dacopnr3s(mfi_br,'in/br.dac',mtype,ix0,jx0,kx0,nx0)
      mfi_bph=76
      call dacopnr3s(mfi_bph,'in/bph.dac',mtype,ix0,jx0,kx0,nx0)
      mfi_bz=77
      call dacopnr3s(mfi_bz,'in/bz.dac',mtype,ix0,jx0,kx0,nx0)
      do n=1,ndi
        read(mfi_t,end=9900) t
        read(mfi_ro) ro
        read(mfi_pr) pr
        read(mfi_vr) vr
        read(mfi_vph) vph
        read(mfi_vz) vz
        read(mfi_br) br
        read(mfi_bph) bph
        read(mfi_bz) bz
      enddo
9900  continue

      endif

c----------------------------------------------------------------------|
c     data output
      mf_r=11
      call dacdef1d(mf_r,'r.dac',6,ix)
      write(mf_r) r
      mf_ph=12
      call dacdef1d(mf_ph,'ph.dac',6,jx)
      write(mf_ph) ph
      mf_z=13
      call dacdef1d(mf_z,'z.dac',6,kx)
      write(mf_z) z
      call dacputparamd(mf_params,'gm',gm)

      nd=nd+1
      write(mf_t) t
      write(mf_ro) ro
      write(mf_pr) pr
      write(mf_vr) vr
      write(mf_vph) vph
      write(mf_vz) vz
      write(mf_br) br
      write(mf_bph) bph
      write(mf_bz) bz
      write(6,913) ns,t,nd
      open(mf_out,file='out.txt',status='old',form='formatted'
     &    ,position='append')
            write(mf_out,913) ns,t,nd
      close(mf_out)

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
         call cfl_m(dt,safety,dtmin,merr,gm,ro,pr,vr,vph,vz,br,bph,bz
     &              ,dr,ix,dph,jx,dz,kx)
         if (merr.ne.0) goto 9999
         tp = t
         t = t+dt

c----------------------------------------------------------------------|
c     solve hydrodynamic equations

         qav=3.d0
c        call mlw_m_c(ro,pr,vr,vph,vz,br,bph,bz,dt,qav,gm
c    &          ,r,rm,dr,drm,ix,dph,dphm,jx,dz,dzm,kx)

         call engine(ro,pr,vr,vph,vz,br,bph,bz,r,rm
     &          ,gm,mu,dt,dr,drm,dph,dphm,dz,dzm,ix,jx,kx,mfdim)

         call bnd(margin,ro,pr,vr,vph,vz,br,bph,bz,ix,jx,kx,mfdim)
c        floor=1.d-9
c        call chkdav(n_floor,ro,vr,vph,vz,floor,ix,jx,kx)
c        call chkdav(n_floor,pr,vr,vph,vz,floor,ix,jx,kx)


c----------------------------------------------------------------------|
c     data output

         mw=0
         nt1=int(tp/dtout)
         nt2=int(t/dtout)
         if (nt1.lt.nt2) mw=1
         if (mw.ne.0) then
            nd=nd+1
           write(mf_t) t
           write(mf_ro) ro
           write(mf_pr) pr
           write(mf_vr) vr
           write(mf_vph) vph
           write(mf_vz) vz
           write(mf_br) br
           write(mf_bph) bph
           write(mf_bz) bz

           write(6,913) ns,t,nd
           open(mf_out,file='out.txt',status='old',form='formatted'
     &         ,position='append')
                 write(mf_out,913) ns,t,nd
           close(mf_out)
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
           write(mf_vr) vr
           write(mf_vph) vph
           write(mf_vz) vz
           write(mf_br) br
           write(mf_bph) bph
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
