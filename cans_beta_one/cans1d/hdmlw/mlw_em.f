c======================================================================|
      subroutine mlw_em(ey,ez,by,bz,dt,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_em
c
c PURPOSE
c    solve vaccuum Maxwell eqs. by modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    ey,ez(ix): [double] Electric field
c    by,bz(ix): [double] Magnetic field
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2005-4-9 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix),dxi(ix),dxim(ix)
      dimension ey(ix),ez(ix),by(ix),bz(ix)
      dimension eyh(ix),ezh(ix),byh(ix),bzh(ix)
      dimension dey(ix),dez(ix),dby(ix),dbz(ix)
      dimension fx(ix)

c----------------------------------------------------------------------|
c     ready
c----------------------------------------------------------------------|
      do i=1,ix
         dxim(i) = 1.0/dxm(i)
         dxi(i) = 1.0/dx(i)
      enddo

c----------------------------------------------------------------------|
c     initialize dro etc.                                   
c----------------------------------------------------------------------|
      do i=1,ix
         dey(i) = 0.0d0
         dez(i) = 0.0d0
         dby(i) = 0.0d0
         dbz(i) = 0.0d0
      enddo

c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  ey ---
      do i=1,ix
         fx(i)=  bz(i)
      enddo
      call mlwhalf(ey ,eyh ,dey,dt,fx,dxi,dxim,ix)

c---  ez ---
      do i=1,ix
         fx(i)= -by(i)
      enddo
      call mlwhalf(ez ,ezh ,dez,dt,fx,dxi,dxim,ix)

c---  by ---
      do i=1,ix
         fx(i)= -ez(i)
      enddo
      call mlwhalf(by ,byh ,dby,dt,fx,dxi,dxim,ix)

c---  bz ---
      do i=1,ix
         fx(i)=  ey(i)
      enddo
      call mlwhalf(bz ,bzh ,dbz,dt,fx,dxi,dxim,ix)

c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  ey  ---
      do i=1,ix-1
         fx(i)=  bzh(i)
      enddo
      call mlwfull(dey ,dt,fx,dxi,ix)

c---  ez  ---
      do i=1,ix-1
         fx(i)= -byh(i)
      enddo
      call mlwfull(dez ,dt,fx,dxi,ix)

c---  by  ---
      do i=1,ix-1
         fx(i)= -ezh(i)
      enddo
      call mlwfull(dby ,dt,fx,dxi,ix)

c---  bz  ---
      do i=1,ix-1
         fx(i)=  eyh(i)
      enddo
      call mlwfull(dbz ,dt,fx,dxi,ix)

c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         ey(i) = ey(i) +dey(i)
         ez(i) = ez(i) +dez(i)
         by(i) = by(i) +dby(i)
         bz(i) = bz(i) +dbz(i)
      enddo

      return
      end
