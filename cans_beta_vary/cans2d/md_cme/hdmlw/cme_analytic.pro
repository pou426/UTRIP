pro cme_analytic,x,y,time,ro,pr,vx,vy,vz,bx,by,bz,az $
   ,zetac_0=zetac_0,rnu_0=rnu_0,eta_0=eta_0,a0_0=a0_0,rlambda_0=rlambda_0 $
   ,rMsun=rMsun,rr0=rr0,rn0_0=rn0_0

;--- Parameters taken from Stone & Norman (1992)
;--- lambda may be wrong in S&N. So I changed it.

if (n_elements(zetac_0) eq 0) then zetac_0=1.104d11  ; (position) of outflow region
if (n_elements(rnu_0  ) eq 0) then rnu_0=2.42d20     ; entropy
if (n_elements(eta_0  ) eq 0) then eta_0=5.24d-8
if (n_elements(a0_0   ) eq 0) then a0_0=1.5d21       ; magnetic flux
if (n_elements(rlambda_0) eq 0) then rlambda_0=5.54d-11  ; (position) of CME

;--- physical constant

      g0= 6.67d-8     ; Gravitational Constant
      rm_0=1.673d-24  ; Mass of proton

;--- normalization

if (n_elements(rMsun) eq 0) then rMsun=2.0d33   ; Mass of the star
if (n_elements(rr0  ) eq 0) then rr0= 1.d11     ; Radius of the star
if (n_elements(rn0_0) eq 0) then rn0_0=1.d8     ; Number density at the surface

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ix=n_elements(x)
jx=n_elements(y)

ro=fltarr(ix,jx)
pr=fltarr(ix,jx)
vx=fltarr(ix,jx)
vy=fltarr(ix,jx)
vz=fltarr(ix,jx)
bx=fltarr(ix,jx)
by=fltarr(ix,jx)
bz=fltarr(ix,jx)
az=fltarr(ix,jx)

      dxm=mkdelta(x)
      dym=mkdelta(y)


      vv0= sqrt(g0*rMsun/rr0)
      tt0= rr0/vv0
      ro0= rm_0*rn0_0
      pr0= ro0*vv0^2
      bb0= sqrt(pr0)


      zetac = zetac_0/rr0
      rnu = rnu_0/(vv0^2/ro0^(1/3.))
      eta = eta_0*tt0^2
      a0 = a0_0/(sqrt(pr0)*rr0^2)
      rlambda = rlambda_0*rr0

      zlam=rlambda*zetac
      zlam2=sqrt(zlam)

      sj1=sqrt(2/!pi/zlam)*(sin(zlam)/zlam-cos(zlam))
      p0= -zlam2*sj1

      phi=sqrt(eta)*time
      rc=zetac*phi
      rs=phi^(7.d0/6.d0)
      d0=exp(-2.d0/3.d0/eta)

      for i=0,ix-1 do begin
        zeta=x(i)/phi
        zlam=rlambda*zeta
        zlam2=sqrt(zlam)
        sj0=sqrt(2/!pi/zlam)*sin(zlam)
        sj1=sqrt(2/!pi/zlam)*(sin(zlam)/zlam-cos(zlam))
      for j=0,jx-1 do begin
        if (x(i) lt rc) then begin
            pr(i,j)=1.d0/(4.d0*rnu^3*x(i)^4)
            ro(i,j)=1.d0/(rnu^3*x(i)^3)
            vx(i,j)=x(i)/time

            bx(i,j)=2.d0*a0/x(i)^2*(p0+zlam2*sj1)*cos(y(j))
            by(i,j)=-a0*rlambda/x(i)/phi*(zlam2*sj0-sj1/zlam2)*sin(y(j))
            bz(i,j)=a0*rlambda/x(i)/phi*(p0+zlam2*sj1)*sin(y(j))
            aa=(p0+zlam2*sj1)*sin(y(j))^2
            ro(i,j)=ro(i,j)+0.50d0/!pi*a0^2/x(i)^3*p0*(4.d0-zlam^2)*aa
            pr(i,j)=pr(i,j)+0.25d0/!pi*a0^2/x(i)^4*p0*(2.d0-zlam^2)*aa
        endif else begin
            bx(i,j)=0.d0
            by(i,j)=0.d0
            bz(i,j)=0.d0
          if (x(i) lt rs) then begin
            pr(i,j)=(7.d0/6.d0)*d0*eta/phi^4*exp(2.d0/3.d0/eta/zeta^9)
            ro(i,j)=7.d0*d0/(phi^3*zeta^8)*exp(2.d0/3.d0/eta/zeta^9)
            vx(i,j)=x(i)/time
          endif else begin
            pr(i,j)=1.d-20
            ro(i,j)=d0/(x(i)^(26.d0/7.d0)) $
                *exp(2.d0/3.d0/eta/x(i)^(9.d0/7.d0))
            vx(i,j)=0.d0
          endelse
        endelse

        vy(i,j)=0.0d0
        vz(i,j)=0.0d0

      endfor
      endfor

;  vector potential

      az(1,0)=0.

      i=0
      for j=1,jx-1 do begin
        az(i,j) = az(i,j-1) $
                + 0.5*(x(i)*sin(y(j-1))*bx(i,j-1) $
                      +x(i)*sin(y(j  ))*bx(i,j  ))*x(i)*dym(j-1)
      endfor
      for i=1,ix-1 do begin
      for j=0,jx-1 do begin
         az(i,j) = az(i-1,j) $
                - 0.5*(x(i-1)*sin(y(j))*by(i-1,j) $
                      +x(i  )*sin(y(j))*by(i  ,j))*dxm(i-1)
      endfor
      endfor

      for i=0,ix-1 do begin
      for j=0,jx-1 do begin
         az(i,j) = az(i,j)/( x(i)*sin(y(j)) )
      endfor
      endfor

end
