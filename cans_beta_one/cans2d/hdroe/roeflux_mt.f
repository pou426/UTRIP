c======================================================================|
      subroutine roeflux_mt(fro,frx,fry,fby,cs2
     &              ,row,vxw,vyw,bxw,byw,ix,jx)
c======================================================================|
c
c NAME  roeflux_mt
c
c PURPOSE
c    derive numerical flux by solving the linearized Riemann problem
c        * isothermal MHD
c
c INPUTS & OUTPUTS
c    None
c
c OUTPUTS
c    fro(ix): [double] density flux
c    frx(ix): [double] momentum flux
c    fry(ix): [double] momentum flux
c    fby(ix): [double] magnetic field flux
c
c INPUTS
c    row(ix,2): [double] density at cell boundary
c    vxw(ix,2): [double] velocity at cell boundary
c    vyw(ix,2): [double] velocity at cell boundary
c    byw(ix,2): [double] magnetic field at cell boundary
c    cs2: [double] square of sound speed
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on N. Fukuda's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension row(ix,jx,2)
      dimension vxw(ix,jx,2),vyw(ix,jx,2)
      dimension bxw(ix,jx,2),byw(ix,jx,2)
      dimension fro(ix,jx),frx(ix,jx),fry(ix,jx),fby(ix,jx)

c----------------------------------------------------------------------|

      pi = acos(-1.0d0)
      pi4=4.0d0*pi
      pi4i=1.0d0/pi4
      pi8i=5.0d-1*pi4i

      do j=1,jx-1
      do i=1,ix-1
         rol=row(i,j,1)
         vxl=vxw(i,j,1)
         vyl=vyw(i,j,1)
         bxl=bxw(i,j,1)
         byl=byw(i,j,1)
         ror=row(i,j,2)
         vxr=vxw(i,j,2)
         vyr=vyw(i,j,2)
         bxr=bxw(i,j,2)
         byr=byw(i,j,2)

         prl=rol*cs2
         prr=ror*cs2

c
c-----barred variables
        sr0=sqrt(rol)
        sr1=sqrt(ror)
        sri=1.0d0/(sr0+sr1)
        robar=sr0*sr1
        vxbar=(sr0*vxl+sr1*vxr)*sri
        vybar=(sr0*vyl+sr1*vyr)*sri
        bxbar=(sr0*bxr+sr1*bxl)*sri
        bybar=(sr0*byr+sr1*byl)*sri
        byave=(byl+byr)/2.0d0
c-----characteristic speed
        cs2x=cs2+((byl-byr)**2)*sri**2*pi8i
        astar2=cs2x+(bxbar**2+bybar**2)*pi4i/robar
        ca2=bxbar**2/(pi4*robar)
        bybz=(bybar**2)*pi4i/robar
        csca=cs2x-ca2
        cfcs=sqrt(csca**2+bybz*(2.0d0*(cs2x+ca2)+bybz))
        cfast2=0.5d0*(astar2+cfcs)
        cslow2=cs2x*ca2/cfast2
        cfast=sqrt(cfast2)
        cslow=sqrt(cslow2)
        ca=sqrt(ca2)
        csx=sqrt(cs2x)
c-----for singular points
        eps=1.0d-12
        sgr=bybar**2-eps
        sp=0.5d0+sign(0.5d0,sgr)
        betay=sp*bybar*sqrt(1.0d0/(bybar**2+1.0d0-sp))
     1         +sqrt(0.5d0)*(1.0d0-sp)
        betaz=  sqrt(0.5d0)*(1.0d0-sp)
        eps2=1.0d-6
        sgr2=cfcs-eps2
        sp2=0.5d0+sign(0.5d0,sgr2)
c  cfca = c_fast ^2 - c_alfven ^2
c  cfa  = c_fast ^2 - c_s ^2
        cfca=(bybz*(2.0*(ca2+cs2x)+bybz))/(abs(csca)+cfcs)
     1      +max(csca,0.0d0)
        cfa=(bybz*(2.0*(ca2+cs2x)+bybz))/(abs(csca)+cfcs)
     1      -min(csca,0.0d0)
        alphf=sp2*sqrt(cfca/(cfcs+1.0d0-sp2))
     1          +1.0d0-sp2
        alphs=sp2*sqrt((cfa)/(cfcs+1.0d0-sp2))

        sgnbx=sign(1.0d0,bxbar)
c-----eigen values
        eeps=(vxr-vxl+abs(vxr-vxl))*2.5d-1 
        elpf=-max(abs(vxbar+cfast),eeps)
        elmf=-max(abs(vxbar-cfast),eeps)
        elps=-max(abs(vxbar+cslow),eeps)
        elms=-max(abs(vxbar-cslow),eeps)
        elpa=-max(abs(vxbar+ca),eeps)
        elma=-max(abs(vxbar-ca),eeps)
        elze=-max(abs(vxbar),eeps)
c-----amplitude;w's
        dro=ror-rol
        du21=robar*(vxr-vxl)
        du31=robar*(vyr-vyl)
        du5=byr-byl
        du6=0.d0
c
        t1=betay*du5+betaz*du6
        t3=betaz*du5-betay*du6

c
        s1=du21
        s2=betay*du31
        s3=betaz*du31
c
        p11=alphs*cfast*sqrt(pi4/robar)
        p12=-alphf*cs2x/cfast*sqrt(pi4/robar)
        p21=alphf
        p22=alphs
c
        q11=alphf*cfast
        q12=alphs*cslow
        q21=-alphs*ca*sgnbx
        q22=alphf*csx*sgnbx
c
        detp=p11*p22-p12*p21
        detq=q11*q22-q12*q21
c
        wpf=0.5d0*((p22*t1-p12*dro)/detp+(q22*s1-q12*s2)/detq)
        wmf=0.5d0*((p22*t1-p12*dro)/detp-(q22*s1-q12*s2)/detq)
        wps=0.5d0*((-p21*t1+p11*dro)/detp+(-q21*s1+q11*s2)/detq)
        wms=0.5d0*((-p21*t1+p11*dro)/detp-(-q21*s1+q11*s2)/detq)
        wpa=0.5d0*(sqrt(robar*pi4i)*t3-sgnbx*s3)
        wma=0.5d0*(sqrt(robar*pi4i)*t3+sgnbx*s3)
        wze=dro-alphf*(wpf+wmf)-alphs*(wps+wms)

c----- flux_l
        fluxlro=rol*vxl
        fluxlrx=rol*vxl*vxl+prl+(-bxbar**2+byl**2)*pi8i
        fluxlry=rol*vxl*vyl-bxbar*byl*pi4i
        fluxlby=vxl*byl-vyl*bxbar
c----- flux_r
        fluxrro=ror*vxr
        fluxrrx=ror*vxr*vxr+prr+(-bxbar**2+byr**2)*pi8i
        fluxrry=ror*vxr*vyr-bxbar*byr*pi4i
        fluxrby=vxr*byr-vyr*bxbar

c-----components of the eigen vectors
        rpfro=alphf
        rpfrx=alphf*(vxbar+cfast)
        rpfry=alphf*vybar-alphs*betay*ca*sgnbx
        rpfby=alphs*betay*cfast*sqrt(pi4/robar)
c
        rmfro=alphf
        rmfrx=alphf*(vxbar-cfast)
        rmfry=alphf*vybar+alphs*betay*ca*sgnbx
        rmfby=rpfby
c
        rpsro=alphs
        rpsrx=alphs*(vxbar+cslow)
        rpsry=alphs*vybar+csx*sgnbx*alphf*betay
        rpsby=-sqrt(pi4/robar)*cs2x*alphf*betay/cfast
c
        rmsro=alphs
        rmsrx=alphs*(vxbar-cslow)
        rmsry=alphs*vybar-csx*sgnbx*alphf*betay
        rmsby=rpsby
c
        rparo=0.0d0
        rparx=0.0d0
        rpary=-sgnbx*betaz
        rpaby=sqrt(pi4/robar)*betaz
c
        rmaro=0.0d0
        rmarx=0.0d0
        rmary=-rpary
        rmaby=rpaby
c
        rzero=1.0d0
        rzerx=vxbar
        rzery=vybar
        rzeby=0.0d0
c
c-----computation of f(i+1/2,j)
        fro(i,j)=0.5d0*(fluxlro+fluxrro
     1            +elpf*wpf*rpfro +elmf*wmf*rmfro
     5            +elpa*wpa*rparo +elma*wma*rmaro
     3            +elps*wps*rpsro +elms*wms*rmsro +elze*wze*rzero)

        frx(i,j)=0.5d0*(fluxlrx+fluxrrx
     1            +elpf*wpf*rpfrx +elmf*wmf*rmfrx
     5            +elpa*wpa*rparx +elma*wma*rmarx
     3            +elps*wps*rpsrx +elms*wms*rmsrx +elze*wze*rzerx)

        fry(i,j)=0.5d0*(fluxlry+fluxrry
     1            +elpf*wpf*rpfry +elmf*wmf*rmfry
     3            +elps*wps*rpsry +elms*wms*rmsry
     5            +elpa*wpa*rpary +elma*wma*rmary +elze*wze*rzery)

        fby(i,j)=0.5d0*(fluxlby+fluxrby
     1            +elpf*wpf*rpfby +elmf*wmf*rmfby
     3            +elps*wps*rpsby +elms*wms*rmsby
     5            +elpa*wpa*rpaby +elma*wma*rmaby +elze*wze*rzeby)

       enddo
       enddo


       return
       end
