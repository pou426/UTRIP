!=====================================================================
subroutine hlld_flux(row,prw,vxw,vyw,vzw,bx,byw,bzw,gm,ix,jx,kx &
     ,fro,fee,frx,fry,frz,fby,fbz)
!=====================================================================
!
! Name :: hlld_flux
!
! Purpose
!
! Reference
!    T.Miyoshi, K.Kusano, JCP, 208, 315 (2005)
!
! Input
!   ix,jx,kx :: array size
!   row,prw,vxw,vyw,vzw,byw,bzw
!            :: primitive value array(ix,jx,kx,2)
!               qqw(i,j,k,1) = qq(i,j,k)'s left state
!               qqw(i,j,k,2) = qq(i,j,k)'s right state
!   bx       :: magnetic field normal to cell surface
!   gm       :: specific heat rate
! Output
!   f**      :: numerical flux
!=====================================================================

  implicit none

!----> INPUT
  integer,intent(in) :: ix,jx,kx ! array size
  real(8),intent(in) :: gm       ! specific heat rate

! primitive variables :: 1 = left state , 2 = right state
! 
  real(8),dimension(ix,jx,kx,2) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2) :: byw,bzw

  real(8),dimension(ix,jx,kx) :: bx ! magnetic field 
                                    ! parallel to x-direction
!----< INPUT

!----> OUTPUT

  real(8),dimension(ix,jx,kx) :: fro,fee,frx,fry,frz ! numerical flux
  real(8),dimension(ix,jx,kx) :: fby,fbz

!----< OUTPUT

!----- U -----
! qql :: left state
! qqr :: right state

  real(8) :: rol,vxl,vyl,vzl,byl,bzl,ptl,eel 
  real(8) :: ror,vxr,vyr,vzr,byr,bzr,ptr,eer 
  
  real(8) :: rxl,ryl,rzl 
  real(8) :: rxr,ryr,rzr
  
  real(8) :: bxs,bxsq
  
  real(8) :: pbl,pbr,ptst,prl,prr
!----- U* ----
! qqlst :: left state
! qqrst :: right state

  real(8) :: rolst,vxlst,vylst,vzlst,bylst,bzlst,eelst
  real(8) :: rorst,vxrst,vyrst,vzrst,byrst,bzrst,eerst
  
  real(8) :: rxlst,rylst,rzlst
  real(8) :: rxrst,ryrst,rzrst
!----- U** ---
! qqldst :: left state
! qqrdst :: right state

  real(8) :: roldst,vxldst,vyldst,vzldst,byldst,bzldst,eeldst
  real(8) :: rordst,vxrdst,vyrdst,vzrdst,byrdst,bzrdst,eerdst
  
  real(8) :: rxldst,ryldst,rzldst
  real(8) :: rxrdst,ryrdst,rzrdst

!----- flux ---
! fqql :: left physical flux
! fqqr :: right physical flux
! fluxqq :: intemediate HLLD flux (OUTPUT)

  real(8) :: frol,frxl,fryl,frzl
  real(8) :: fbyl,fbzl,feel
  
  real(8) :: fror,frxr,fryr,frzr
  real(8) :: fbyr,fbzr,feer

!----- wave speed ---
! sl :: left-going fastest signal velocity
! sr :: right-going fastest signal velocity
! sm :: contact discontinuity velocity
! slst :: left-going alfven velocity
! srst :: right-going alfven velocity
  real(8) :: sm,sl,sr,slst,srst

! cfl :: left-state Fast wave velocity
! cfr :: right-sate Fast wave velocity
! csl :: left-state Slow wave velocity
! csr :: right-state Slow wave velocity
! cal :: left-state Alfven wave velocity
! car :: right-state Alfven wave velocity
! cacl :: left-state Sound wave velocity
! cacr :: right-state Sound wave velocity
  real(8) :: cfl,cfr,csl,csr,cal,car,cacl,cacr 

!--------------------
! temporary variables

  real(8) :: gmpl,gmpr
  real(8) :: gpbl,gpbr
  real(8) :: cfmax,maxspd
  
  real(8) :: sdl,sdr,sdml,sdmr
  real(8) :: vdbstl,vdbstr
  real(8) :: sqrtrol,sqrtror
  real(8) :: invsumro
  real(8) :: signbx,inverse_sqrtro,msqrtro
  
  real(8) :: bsq,sqrt_bsq
  real(8) :: bsqr,bsql
  real(8) :: temp

  real(8) :: signBy,signBz
  real(8) :: temp_fst,sqrt_fst,temp_fa,sqrt_fa
  real(8) :: alpha_1,alpha_2,alpha_3  
  real(8) :: beta_y,beta_z
  real(8) :: mu,pi
  
  integer :: i,j,k  

      pi = acos(-1.0d0)
      mu= 4.d0*pi
  

  do k=1,kx-1
     do j=1,jx-1
        do i=1,ix-1

!----- Step 0. ----------------------------------------------------------|
! set L/R-state
!
           bxs = bx(i,j,k)
           bxsq = bxs**2
!---- Left state

           rol = row(i,j,k,1)
           vxl = vxw(i,j,k,1)
           vyl = vyw(i,j,k,1)
           vzl = vzw(i,j,k,1)
           byl = byw(i,j,k,1)
           bzl = bzw(i,j,k,1)
           
           bsql = bxs**2+byl**2+bzl**2
!          pbl = 0.5d0*(bxs**2 + byl**2 + bzl**2)
           pbl = (bxs**2 + byl**2 + bzl**2)/(2.d0*mu)
           prl = prw(i,j,k,1)
           ptl = prl + pbl
           
           eel = prl/(gm-1.0d0) &
                + 0.5d0*rol*(vxl**2+vyl**2+vzl**2) &
                + pbl
           rxl = rol*vxl
           ryl = rol*vyl
           rzl = rol*vzl
!---- Right state

           ror = row(i,j,k,2)
           vxr = vxw(i,j,k,2)
           vyr = vyw(i,j,k,2)
           vzr = vzw(i,j,k,2)
           byr = byw(i,j,k,2)
           bzr = bzw(i,j,k,2)
           
           bsqr = bxs**2+byr**2+bzr**2
!          pbr = (bxs**2 + byr**2 + bzr**2)
           pbr = (bxs**2 + byr**2 + bzr**2)/(2.d0*mu)
           prr = prw(i,j,k,2)
           
           ptr = prr + pbr
           eer = prr/(gm-1.0d0) &
                + 0.5d0*(ror*(vxr**2+vyr**2+vzr**2)) &
                + pbr
           
           rxr = ror*vxr
           ryr = ror*vyr
           rzr = ror*vzr

!----- Step 1. ----------------------------------------------------------|
! Compute wave left & right wave speed
!
           
           gmpl = gm*prl
           gmpr = gm*prr
           
           gpbl = gmpl+2.0d0*pbl
           gpbr = gmpr+2.0d0*pbr
           
           cfl = sqrt((gpbl + sqrt(gpbl**2-4.0d0*gmpl*bxsq/mu))/(2.0d0*rol))
           cfr = sqrt((gpbr + sqrt(gpbr**2-4.0d0*gmpr*bxsq/mu))/(2.0d0*ror))
           
           cfmax = max(cfl,cfr)
           
           if ( vxl <= vxr) then
              sl = vxl - cfmax
              sr = vxr + cfmax
           else
              sl = vxr - cfmax
              sr = vxl + cfmax
           endif
           
           maxspd = max(dabs(sl),dabs(sr))

!----- Step 2. ----------------------------------------------------------|
! compute L/R fluxs
!

! Left value
           frol = rxl
           frxl = rxl*vxl + ptl -bxsq/mu
           fryl = rol*vxl*vyl - bxs*byl/mu
           frzl = rol*vxl*vzl - bxs*bzl/mu
           feel = vxl*(eel + ptl -bxsq/mu) - bxs*(vyl*byl + vzl*bzl)/mu
           fbyl = byl*vxl - bxs*vyl
           fbzl = bzl*vxl - bxs*vzl

! Right value
           fror = rxr
           frxr = rxr*vxr + ptr -bxsq/mu
           fryr = ror*vxr*vyr - bxs*byr/mu
           frzr = ror*vxr*vzr - bxs*bzr/mu
           feer = vxr*(eer + ptr -bxsq/mu) - bxs*(vyr*byr + vzr*bzr)/mu
           fbyr = byr*vxr - bxs*vyr
           fbzr = bzr*vxr - bxs*vzr

!----- Step 3. ----------------------------------------------------------|
! return upwind flux
!

           if (sl >= 0.0d0) then
              fro(i,j,k) = frol
              frx(i,j,k) = frxl
              fry(i,j,k) = fryl
              frz(i,j,k) = frzl
              fee(i,j,k) = feel
              fby(i,j,k) = fbyl
              fbz(i,j,k) = fbzl
              cycle
           endif
           
           if (sr <= 0.0d0) then
              fro(i,j,k) = fror
              frx(i,j,k) = frxr
              fry(i,j,k) = fryr
              frz(i,j,k) = frzr
              fee(i,j,k) = feer
              fby(i,j,k) = fbyr
              fbz(i,j,k) = fbzr
              cycle
           endif

!----- Step 4. ----------------------------------------------------------|
! compute middle and alfven wave
!
           sdl = sl - vxl
           sdr = sr - vxr
           
           sm = (sdr*ror*vxr - sdl*rol*vxl - ptr + ptl) &
                /(sdr*ror - sdl*rol)
           
           sdml = sl - sm
           sdmr = sr - sm
           
!----- Step 5. ----------------------------------------------------------|
! compute intermediate states
!

!              ptst = ptl + rol*sdl*(sdl-sdml)
           ptst = (sdr*ror*ptl-sdl*rol*ptr+rol*ror*sdr*sdl*(vxr-vxl)) &
                /(sdr*ror-sdl*rol)

!----- Step 5A. ----------------------------------------------------------|
! compute Ul*
!

           temp_fst = dabs(rol*sdl*sdml - bxsq/mu)
           sqrt_fst = sqrt(temp_fst)
           temp_fa = dabs(rol*sdl**2 - bxsq/mu)
           sqrt_fa = sqrt(temp_fa)

           bsq = byl**2+bzl**2
           sqrt_bsq = sqrt(bsq/mu)
           
           if(temp_fst .eq. 0.0d0)then
              alpha_1 = 0.0d0
              alpha_2 = 0.0d0
              alpha_3 = 0.0d0
           else
              alpha_1 = sqrt_bsq/sqrt_fst
              alpha_2 = (sm-vxl)/sqrt_fst
              alpha_3 = sqrt_fa/sqrt_fst
           endif

           signBy = sign(1.0d0,byl)
           signBz = sign(1.0d0,bzl)
           if(sqrt_bsq .eq. 0.0d0)then
              beta_z = 0.0d0
              beta_y = 0.0d0
           else
              beta_y = byl/sqrt(mu)/sqrt_bsq
              beta_z = bzl/sqrt(mu)/sqrt_bsq
           endif

           if(sqrt_fst .eq. 0.0d0)then
              rolst = rol

              vxlst = vxl
              rxlst = rolst*vxlst 

              vylst = vyl
              vzlst = vzl
              rylst = rolst*vylst
              rzlst = rolst*vzlst
              
              bylst = byl
              bzlst = bzl
              eelst = eel
           else
              rolst = rol*sdl/sdml

              rxlst = rolst*sm
              vxlst = sm

              vylst = vyl-bxs*beta_y*alpha_1*alpha_2/rol/mu
              vzlst = vzl-bxs*beta_z*alpha_1*alpha_2/rol/mu
              rylst = rolst*vylst
              rzlst = rolst*vzlst
              bylst = beta_y*alpha_1*alpha_3*sqrt_fa*sqrt(mu)
              bzlst = beta_z*alpha_1*alpha_3*sqrt_fa*sqrt(mu)
           
              vdbstl = (rxlst*bxs+rylst*bylst+rzlst*bzlst)/rolst
              eelst = (sdl*eel - ptl*vxl + ptst*sm + bxs*(vxl*bxs+vyl*byl+vzl*bzl-vdbstl)/mu)/sdml
           endif
!----- Step 5B. ----------------------------------------------------------|
! compute Ur*
!
           temp_fst = dabs(ror*sdr*sdmr - bxsq/mu)
           sqrt_fst = sqrt(temp_fst)
           temp_fa = dabs(ror*sdr**2 - bxsq/mu)
           sqrt_fa = sqrt(temp_fa)
           
           bsq = byr**2+bzr**2
           sqrt_bsq = sqrt(bsq/mu)
           
           if(temp_fst .eq. 0.0d0)then
              alpha_1 = 0.0d0
              alpha_2 = 0.0d0
              alpha_3 = 0.0d0
           else
              alpha_1 = sqrt_bsq/sqrt_fst
              alpha_2 = (sm-vxr)/sqrt_fst
              alpha_3 = sqrt_fa/sqrt_fst
           endif

           signBy = sign(1.0d0,byr)
           signBz = sign(1.0d0,bzr)
           if(sqrt_bsq .eq. 0.0d0)then
              beta_z = 0.0d0
              beta_y = 0.0d0
           else
              beta_y = byr/sqrt(mu)/sqrt_bsq
              beta_z = bzr/sqrt(mu)/sqrt_bsq
           endif

           if(sqrt_fst .eq. 0.0d0)then
              rorst = ror

              vxrst = vxr
              rxrst = rorst*vxrst 

              vyrst = vyr
              vzrst = vzr
              ryrst = rorst*vyrst
              rzrst = rorst*vzrst
              
              byrst = byr
              bzrst = bzr
              eerst = eer
           else
              rorst = ror*sdr/sdmr

              rxrst = rorst * sm
              vxrst = sm

              vyrst = vyr-bxs*beta_y*alpha_1*alpha_2/ror/mu
              vzrst = vzr-bxs*beta_z*alpha_1*alpha_2/ror/mu
              ryrst = rorst*vyrst
              rzrst = rorst*vzrst

              byrst = beta_y*alpha_1*alpha_3*sqrt_fa*sqrt(mu)
              bzrst = beta_z*alpha_1*alpha_3*sqrt_fa*sqrt(mu)

              vdbstr = (rxrst*bxs+ryrst*byrst+rzrst*bzrst)/rorst
              eerst = (sdr*eer - ptr*vxr + ptst*sm + bxs*(vxr*bxs+vyr*byr+vzr*bzr-vdbstr)/mu)/sdmr
           endif
!----- Step 5C. ----------------------------------------------------------|
! compute Ul** and Ur**
!
           sqrtrol = sqrt(rolst)
           sqrtror = sqrt(rorst)
           
           slst = sm - dabs(bxs)/sqrtrol/sqrt(mu)
           srst = sm + dabs(bxs)/sqrtror/sqrt(mu)

           if (bxs .eq. 0.0d0)then
              roldst = rolst
              rxldst = rxlst
              ryldst = rylst
              rzldst = rzlst
              vxldst = vxlst
              vyldst = vylst
              vzldst = vzlst
              eeldst = eelst
              byldst = bylst
              bzldst = bzlst
              
              rordst = rorst
              rxrdst = rxrst
              ryrdst = ryrst
              rzrdst = rzrst
              vxrdst = vxrst
              vyrdst = vyrst
              vzrdst = vzrst
              eerdst = eerst
              byrdst = byrst
              bzrdst = bzrst
           else
              invsumro = 1.0d0/(sqrtrol + sqrtror)
              signbx = sign(1.0d0,bxs)
              
              roldst = rolst
              rordst = rorst
              
              rxldst = rxlst
              rxrdst = rxrst
              vxldst = vxlst
              vxrdst = vxrst
              
              temp = invsumro*(sqrtrol*vylst + sqrtror*vyrst + signbx*(byrst-bylst)/mu)
              vyldst = temp
              vyrdst = temp
              ryldst = roldst * temp
              ryrdst = rordst * temp
              
              temp = invsumro*(sqrtrol*vzlst + sqrtror*vzrst + signbx*(bzrst-bzlst)/mu)
              vzldst = temp
              vzrdst = temp
              rzldst = roldst * temp
              rzrdst = rordst * temp
              
              temp = invsumro*(sqrtrol*byrst + sqrtror*bylst  &
                   + signbx*sqrtrol*sqrtror*(vyrst-vylst))
              byldst = temp
              byrdst = temp
              
              temp = invsumro*(sqrtrol*bzrst + sqrtror*bzlst &
                   + signbx*sqrtrol*sqrtror*(vzrst-vzlst))
              
              bzldst = temp
              bzrdst = temp
              
              temp = sm*bxs + vyldst*byldst + vzldst*bzldst
              eeldst = eelst - sqrtrol*signbx*(vdbstl - temp)
              eerdst = eerst + sqrtror*signbx*(vdbstr - temp)
           endif

!----- Step 6. ----------------------------------------------------------|
! compute flux
!

           if (slst >= 0.0d0)then
              fro(i,j,k) = frol + sl*(rolst - rol)
              frx(i,j,k) = frxl + sl*(rxlst - rxl)
              fry(i,j,k) = fryl + sl*(rylst - ryl)
              frz(i,j,k) = frzl + sl*(rzlst - rzl)
              fee(i,j,k) = feel + sl*(eelst - eel)
              fby(i,j,k) = fbyl + sl*(bylst - byl)
              fbz(i,j,k) = fbzl + sl*(bzlst - bzl)
           else if(sm >= 0.0d0) then
              temp = slst - sl
              fro(i,j,k) = frol - sl*rol - temp*rolst + slst*roldst
              frx(i,j,k) = frxl - sl*rxl - temp*rxlst + slst*rxldst
              fry(i,j,k) = fryl - sl*ryl - temp*rylst + slst*ryldst
              frz(i,j,k) = frzl - sl*rzl - temp*rzlst + slst*rzldst
              fee(i,j,k) = feel - sl*eel - temp*eelst + slst*eeldst
              fby(i,j,k) = fbyl - sl*byl - temp*bylst + slst*byldst
              fbz(i,j,k) = fbzl - sl*bzl - temp*bzlst + slst*bzldst
           else if (srst > 0.0d0) then
              temp = srst - sr
              fro(i,j,k) = fror - sr*ror - temp*rorst + srst*rordst
              frx(i,j,k) = frxr - sr*rxr - temp*rxrst + srst*rxrdst
              fry(i,j,k) = fryr - sr*ryr - temp*ryrst + srst*ryrdst
              frz(i,j,k) = frzr - sr*rzr - temp*rzrst + srst*rzrdst
              fee(i,j,k) = feer - sr*eer - temp*eerst + srst*eerdst
              fby(i,j,k) = fbyr - sr*byr - temp*byrst + srst*byrdst
              fbz(i,j,k) = fbzr - sr*bzr - temp*bzrst + srst*bzrdst
           else
              fro(i,j,k) = fror + sr*(rorst - ror)
              frx(i,j,k) = frxr + sr*(rxrst - rxr)
              fry(i,j,k) = fryr + sr*(ryrst - ryr)
              frz(i,j,k) = frzr + sr*(rzrst - rzr)
              fee(i,j,k) = feer + sr*(eerst - eer)
              fby(i,j,k) = fbyr + sr*(byrst - byr)
              fbz(i,j,k) = fbzr + sr*(bzrst - bzr)
           endif
           
        end do
     end do
  enddo
end subroutine hlld_flux
