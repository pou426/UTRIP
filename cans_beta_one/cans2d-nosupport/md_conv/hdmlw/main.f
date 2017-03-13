c======================================================================|
c     array definitions
c======================================================================|
      implicit real*8 (a-h,o-z)
      parameter (ix=145,jx=64)

      dimension x(ix),xm(ix)
      dimension dx(ix),dxm(ix)
      dimension y(jx),ym(jx)
      dimension dy(jx),dym(jx)

      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx)
      dimension gx(ix,jx),gy(ix,jx)
      dimension gxm(ix,jx),gym(ix,jx)
      dimension tem0(jx),den0(jx)
      dimension rkap(ix,jx),rkapm(ix,jx),visc(ix,jx),viscm(ix,jx)
      dimension hcs(ix,jx),hcsm(ix,jx)



c======================================================================|
c     prologue
c======================================================================|

c----------------------------------------------------------------------|
c   parameters

      margin=2

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

       call dacputparamc(mf_params,'comment','cans2d md_cv')
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

      tend=100.
      dtout=1.
      nstop=3500000
c     dtout=1.0d-7
c     nstop=6

c----------------------------------------------------------------------|
c   setup numerical model (grid, initial conditions, etc.)

      call model(ro,pr,vx,vy,gm
     &           ,gx,gxm,gy,gym,tem0,den0
     &     ,margin,x,ix,y,jx,rkap,rkapm,visc,viscm,hcs,hcsm,gasr)

      call pertub(vy,x,ix,y,jx,margin)

      teb1=tem0(3)
      tebx=tem0(jx-2)

      call grdrdy(dy,ym,dym,y,jx)

      call bnd(margin,ro,pr,vx,vy,ix,jx,gasr,teb1,tebx)

      floor=1.d-9
      call chkdav(n_floor,ro,vx,vy,floor,ix,jx)
      call chkdav(n_floor,pr,vx,vy,floor,ix,jx)


c----------------------------------------------------------------------|
c     ready

      call grdrdy(dx,xm,dxm,x,ix)
      call grdrdy(dy,ym,dym,y,jx)



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

         call cfl_h_kv(dt,merr,gm,ro,pr,vx,vy
     &           ,rkap,visc,dx,ix,dy,jx)

         if (merr.ne.0) goto 9999
         timep = time
         time = time+dt

c----------------------------------------------------------------------|
c     solve hydrodynamic equations


c                                                      hdmlw - start >>>

         call mlw_h_tfix(ro,pr,vx,vy,dt,gm,gx,gxm,gy,gym
     &           ,x,dx,dxm,ix,y,dy,dym,jx,rkap,rkapm,visc,viscm
     &           ,hcs,hcsm,gasr,teb1,tebx)

c                                                      hdmlw - end   <<<

         call bnd(margin,ro,pr,vx,vy,ix,jx,gasr,teb1,tebx)

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
            write(mf_vy) vy

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

