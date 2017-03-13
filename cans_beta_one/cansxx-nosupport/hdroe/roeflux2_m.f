c----------------------------------------------------------------------|
c     version 1.1 ( 2005/02/08 Yuji SATO )
c   2006. 7. 18. bug-fixed by T. Yokoyama based on T. Hanayama's suggestion
c----------------------------------------------------------------------|

      subroutine roeflux2_m(fro,fee,frx,fry,frz,fbx,fby,fbz,gm,
     &                     row,prw,vxw,vyw,vzw,bxw,byw,bzw,ix,jx,kx)

      implicit double precision (a-h,o-z)

      dimension row(ix,jx,kx,2),prw(ix,jx,kx,2)
      dimension vxw(ix,jx,kx,2),vyw(ix,jx,kx,2),vzw(ix,jx,kx,2)
      dimension bxw(ix,jx,kx,2),byw(ix,jx,kx,2),bzw(ix,jx,kx,2)

      dimension fro(ix,jx,kx),fee(ix,jx,kx)
      dimension frx(ix,jx,kx),fry(ix,jx,kx),frz(ix,jx,kx)
      dimension fbx(ix,jx,kx),fby(ix,jx,kx),fbz(ix,jx,kx)

c----------------------------------------------------------------------|

      pi   = acos(-1.0d0)
      pi4  = 4.0d0*pi
      pi4i = 1.0d0/pi4
      pi8i = 5.0d-1*pi4i

      do k=1,kx
      do j=1,jx
      do i=1,ix
         rhol = row(i,j,k,1)
         vxl  = vxw(i,j,k,1)
         vyl  = vyw(i,j,k,1)
         vzl  = vzw(i,j,k,1)
         bxl  = bxw(i,j,k,1)
         byl  = byw(i,j,k,1)
         bzl  = bzw(i,j,k,1)
         prl  = prw(i,j,k,1)
         rhor = row(i,j,k,2)
         vxr  = vxw(i,j,k,2)
         vyr  = vyw(i,j,k,2)
         vzr  = vzw(i,j,k,2)
         bxr  = bxw(i,j,k,2)
         byr  = byw(i,j,k,2)
         bzr  = bzw(i,j,k,2)
         prr  = prw(i,j,k,2)

c-----roe's variable
      sr0=sqrt(rhol)
      sr1=sqrt(rhor)
      sri=1.0d0/(sr0+sr1)
      rhobar=sr0*sr1
      vxbar=(sr0*vxl+sr1*vxr)*sri
      vybar=(sr0*vyl+sr1*vyr)*sri
      vzbar=(sr0*vzl+sr1*vzr)*sri
      bxbar=(sr0*bxr+sr1*bxl)*sri
      bybar=(sr0*byr+sr1*byl)*sri
      bzbar=(sr0*bzr+sr1*bzl)*sri
      hl=0.5d0*(vxl**2+vyl**2+vzl**2)+gm*prl/((gm-1.0d0)*rhol)
     1  +(bxbar**2+byl**2+bzl**2)/(pi4*rhol)
      hr=0.5d0*(vxr**2+vyr**2+vzr**2)+gm*prr/((gm-1.0d0)*rhor)
     1  +(bxbar**2+byr**2+bzr**2)/(pi4*rhor)
      hbar=(sr0*hl+sr1*hr)*sri

c-----average 
      byave=(byl+byr)/2.0d0
      bzave=(bzl+bzr)/2.0d0

c-----characteristic speed
      delb2=(gm-2.0d0)/(gm-1.0d0)
     1     *((byr-byl)**2+(bzr-bzl)**2)*sri**2*pi8i
      cs2=(gm-1.0d0)*(hbar-0.5d0*(vxbar**2+vybar**2+vzbar**2)
     1   -delb2-(bxbar**2+bybar**2+bzbar**2)*pi4i/rhobar)
      astar2=(gm-1.0d0)
     1      *(hbar-0.5d0*(vxbar**2+vybar**2+vzbar**2)-delb2)
     2      -(gm-2.0d0)*(bxbar**2+bybar**2+bzbar**2)*pi4i/rhobar
      ca2=bxbar**2/(pi4*rhobar)
      cbr2=(bybar**2+bzbar**2)*pi4i/rhobar
      cfast2=0.5d0*(astar2+sqrt(cbr2*(astar2+cs2+ca2)+(cs2-ca2)**2))
      cslow2=cs2*ca2/cfast2
      cfast=sqrt(cfast2)
      cslow=sqrt(cslow2)
      ca=sqrt(ca2)
      cs=sqrt(cs2)

c----- for singular points
      epsi=1.0d-12
      sgr=bybar**2+bzbar**2-epsi
      sp=0.5d0+sign(0.5d0,sgr)
      betay=sp*bybar*sqrt(1.0d0/(bybar**2+bzbar**2+1.0d0-sp))
     1     +sqrt(0.5d0)*(1.0d0-sp)
      betaz=sp*bzbar*sqrt(1.0d0/(bybar**2+bzbar**2+1.0d0-sp))
     1     +sqrt(0.5d0)*(1.0d0-sp)
      eps2=1.0d-12
      sgr2=(bybar**2+bzbar**2)/(pi4*rhobar)+abs(ca2-cs2)-eps2
      sp2=0.5d0+sign(0.5d0,sgr2)
      cfca=max(0.0d0,cfast2-ca2)
      cfcs=max(0.0d0,cfast2-cslow2)
      cfa=max(0.0d0,cfast2-cs2)
      alphf=sp2*sqrt(cfca/(cfcs+1.0d0-sp2))+1.0d0-sp2
      alphs=sp2*sqrt(cfa/(cfcs+1.0d0-sp2))+1.0d0-sp2
      sgnbx=sign(1.0d0,bxbar)

c----- eigen value & entropy condition
      eeps=(vxr-vxl+abs(vxr-vxl))*2.5d-1 
      elpf=-max(abs(vxbar+cfast),eeps)
      elmf=-max(abs(vxbar-cfast),eeps)
      elps=-max(abs(vxbar+cslow),eeps)
      elms=-max(abs(vxbar-cslow),eeps)
      elpa=-max(abs(vxbar+ca),eeps)
      elma=-max(abs(vxbar-ca),eeps)
      elze=-max(abs(vxbar),eeps)

c----- amplitude;w's
      drho=rhor-rhol
      du21=rhobar*(vxr-vxl)
      du31=rhobar*(vyr-vyl)
      du41=rhobar*(vzr-vzl)
      du6=byr-byl
      du5=bzr-bzl
      t1=betay*du6+betaz*du5
      t2=(prr-prl+(byave*du6+bzave*du5)*pi4i
     1  +(gm-2.0d0)*(bybar*du6+bzbar*du5)*pi4i)/(gm-1.0d0)
      t3=betaz*du6-betay*du5
      s1=du21
      s2=betay*du31+betaz*du41
      s3=betaz*du31-betay*du41
      p11=alphs*cfast*sqrt(pi4/rhobar)
      p12=-alphf*cs2/cfast*sqrt(pi4/rhobar)
      p21=alphf*(cfast2-cs2*(gm-2.0d0)/(gm-1.0d0))
      p22=alphs*(cslow2-cs2*(gm-2.0d0)/(gm-1.0d0))
      q11=alphf*cfast
      q12=alphs*cslow
      q21=-alphs*ca*sgnbx
      q22=alphf*cs*sgnbx
      detp=p11*p22-p12*p21
      detq=q11*q22-q12*q21
      wpf=0.5d0*((p22*t1-p12*t2)/detp+(q22*s1-q12*s2)/detq)
      wmf=0.5d0*((p22*t1-p12*t2)/detp-(q22*s1-q12*s2)/detq)
      wps=0.5d0*((-p21*t1+p11*t2)/detp+(-q21*s1+q11*s2)/detq)
      wms=0.5d0*((-p21*t1+p11*t2)/detp-(-q21*s1+q11*s2)/detq)
      wpa=0.5d0*(sqrt(rhobar*pi4i)*t3-sgnbx*s3)
      wma=0.5d0*(sqrt(rhobar*pi4i)*t3+sgnbx*s3)
      wze=drho-alphf*(wpf+wmf)-alphs*(wps+wms)

c----- flux
      fluxlro=rhol*vxl
      fluxlrx=rhol*vxl*vxl+prl+(-bxbar**2+byl**2+bzl**2)*pi8i
      fluxlry=rhol*vxl*vyl-bxbar*byl*pi4i
      fluxlrz=rhol*vxl*vzl-bxbar*bzl*pi4i
      fluxlby=vxl*byl-vyl*bxbar
      fluxlbz=vxl*bzl-vzl*bxbar
      fluxlee=rhol*vxl*hl-bxbar*(bxbar*vxl+byl*vyl+bzl*vzl)*pi4i
      fluxrro=rhor*vxr
      fluxrrx=rhor*vxr*vxr+prr+(-bxbar**2+byr**2+bzr**2)*pi8i
      fluxrry=rhor*vxr*vyr-bxbar*byr*pi4i
      fluxrrz=rhor*vxr*vzr-bxbar*bzr*pi4i
      fluxrby=vxr*byr-vyr*bxbar
      fluxrbz=vxr*bzr-vzr*bxbar
      fluxree=rhor*vxr*hr-bxbar*(bxbar*vxr+byr*vyr+bzr*vzr)*pi4i

c----- components of the eigen vectors
      rpfro=alphf
      rpfrx=alphf*(vxbar+cfast)
      rpfry=alphf*vybar-alphs*betay*ca*sgnbx
      rpfrz=alphf*vzbar-alphs*betaz*ca*sgnbx
      rpfby=alphs*betay*cfast*sqrt(pi4/rhobar)
      rpfbz=alphs*betaz*cfast*sqrt(pi4/rhobar)
      rpfee=alphf*(0.5d0*(vxbar**2+vybar**2+vzbar**2)
     1           +delb2+cfast*vxbar+cfast2/(gm-1.0d0)
     2           +(cfast2-cs2)*(gm-2.0d0)/(gm-1.0d0))
     3    -alphs*ca*(betay*vybar+betaz*vzbar)*sgnbx

      rmfro=alphf
      rmfrx=alphf*(vxbar-cfast)
      rmfry=alphf*vybar+alphs*betay*ca*sgnbx
      rmfrz=alphf*vzbar+alphs*betaz*ca*sgnbx
      rmfby=rpfby
      rmfbz=rpfbz
      rmfee=alphf*(0.5d0*(vxbar**2+vybar**2+vzbar**2)
     1           +delb2-cfast*vxbar +cfast2/(gm-1.0d0)
     2           +(cfast2-cs2)*(gm-2.0d0)/(gm-1.0d0))
     3    +alphs*ca*(betay*vybar+betaz*vzbar)*sgnbx

      rpsro=alphs
      rpsrx=alphs*(vxbar+cslow)
      rpsry=alphs*vybar+cs*sgnbx*alphf*betay
      rpsrz=alphs*vzbar+cs*sgnbx*alphf*betaz
      rpsby=-sqrt(pi4/rhobar)*cs2*alphf*betay/cfast
      rpsbz=-sqrt(pi4/rhobar)*cs2*alphf*betaz/cfast
      rpsee=alphs*(0.5d0*(vxbar**2+vybar**2+vzbar**2)
     1           +delb2+cslow*vxbar+cslow2/(gm-1.0d0)
     2           +(cslow2-cs2)*(gm-2.0d0)/(gm-1.0d0))
     3    +alphf*cs*(betay*vybar+betaz*vzbar)*sgnbx

      rmsro=alphs
      rmsrx=alphs*(vxbar-cslow)
      rmsry=alphs*vybar-cs*sgnbx*alphf*betay
      rmsrz=alphs*vzbar-cs*sgnbx*alphf*betaz
      rmsby=rpsby
      rmsbz=rpsbz
      rmsee=alphs*(0.5d0*(vxbar**2+vybar**2+vzbar**2)
     1           +delb2-cslow*vxbar+cslow2/(gm-1.0d0)
     2           +(cslow2-cs2)*(gm-2.0d0)/(gm-1.0d0))
     3    -alphf*cs*(betay*vybar+betaz*vzbar)*sgnbx

      rparo=0.0d0
      rparx=0.0d0
      rpary=-sgnbx*betaz
      rparz=sgnbx*betay
      rpaby=sqrt(pi4/rhobar)*betaz
      rpabz=-sqrt(pi4/rhobar)*betay
      rpaee=-(betaz*vybar-betay*vzbar)*sgnbx

      rmaro=0.0d0
      rmarx=0.0d0
      rmary=-rpary
      rmarz=-rparz
      rmaby=rpaby
      rmabz=rpabz
      rmaee=-rpaee

      rzero=1.0d0
      rzerx=vxbar
      rzery=vybar
      rzerz=vzbar
      rzeby=0.0d0
      rzebz=0.0d0
      rzeee=0.5d0*(vxbar**2+vybar**2+vzbar**2)+delb2

c-----computation of f(i+1/2,j)
      fro(i,j,k)=0.5d0*(fluxlro+fluxrro
     1       +elpf*wpf*rpfro +elmf*wmf*rmfro 
     &       +elps*wps*rpsro +elms*wms*rmsro
     2       +elpa*wpa*rparo +elma*wma*rmaro +elze*wze*rzero)

      fee(i,j,k)=0.5d0*(fluxlee+fluxree
     1       +elpf*wpf*rpfee +elmf*wmf*rmfee 
     &       +elps*wps*rpsee +elms*wms*rmsee
     2       +elpa*wpa*rpaee +elma*wma*rmaee +elze*wze*rzeee)

      frx(i,j,k)=0.5d0*(fluxlrx+fluxrrx
     1       +elpf*wpf*rpfrx +elmf*wmf*rmfrx 
     &       +elps*wps*rpsrx +elms*wms*rmsrx
     2       +elpa*wpa*rparx +elma*wma*rmarx +elze*wze*rzerx)

      fry(i,j,k)=0.5d0*(fluxlry+fluxrry
     1       +elpf*wpf*rpfry +elmf*wmf*rmfry 
     &       +elps*wps*rpsry +elms*wms*rmsry
     2       +elpa*wpa*rpary +elma*wma*rmary +elze*wze*rzery)

      frz(i,j,k)=0.5d0*(fluxlrz+fluxrrz
     1       +elpf*wpf*rpfrz +elmf*wmf*rmfrz 
     &       +elps*wps*rpsrz +elms*wms*rmsrz
     2       +elpa*wpa*rparz +elma*wma*rmarz +elze*wze*rzerz)

      fby(i,j,k)=0.5d0*(fluxlby+fluxrby
     1       +elpf*wpf*rpfby +elmf*wmf*rmfby 
     &       +elps*wps*rpsby +elms*wms*rmsby
     2       +elpa*wpa*rpaby +elma*wma*rmaby +elze*wze*rzeby)

      fbz(i,j,k)=0.5d0*(fluxlbz+fluxrbz
     1       +elpf*wpf*rpfbz +elmf*wmf*rmfbz 
     &       +elps*wps*rpsbz +elms*wms*rmsbz
     2       +elpa*wpa*rpabz +elma*wma*rmabz +elze*wze*rzebz)

      fbx(i,j,k)=0.d0
      enddo
      enddo
      enddo

      return
      end

