function shktb_analytic_func,x

common shktb_analytic_param,gm,ro4,pr4,ro1,pr1

cs1=sqrt(gm*pr1/ro1)
cs4=sqrt(gm*pr4/ro4)
cs14=cs1/cs4

power=2.d0*gm/(gm-1.d0)

; Based on "Sawada's Note" eq.(8.274)

y = ( 1.d0- ( (gm-1.d0)/2.d0 *cs14 *(x-1.d0) $
   *sqrt( 2.d0/gm/((gm+1.d0)*x+(gm-1.d0)) ) ) )^power $
   - pr1/pr4*x

return,y
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro shktb_analytic,ro,pr,vx,t,x0,x,gm,ro1,pr1 $
  ,ro4=ro4,pr4=pr4,tolf=tolf,pr21_init=pr210,quiet=quiet

common shktb_analytic_param,gmc,ro4c,pr4c,ro1c,pr1c

if (n_elements(ro4) eq 0) then ro4=1.d0
if (n_elements(pr4) eq 0) then pr4=1.d0
if (n_elements(tolf) eq 0) then tolf=1.d-8 ; Newton-method tolerance
if (n_elements(pr210) eq 0) then pr210=2.d0 ; initial guess

u1=0.d0 & u4=0.d0
gmc=gm & ro4c=ro4 & pr4c=pr4 & ro1c=ro1 & pr1c=pr1

cs1=sqrt(gm*pr1/ro1)
cs4=sqrt(gm*pr4/ro4)

ix=n_elements(x)
ro=fltarr(ix)
pr=fltarr(ix)
vx=fltarr(ix)

;
; shock

pr21 = newton(pr210, 'shktb_analytic_func',tolf=tolf)
if not keyword_set(quiet) then print,'Pr2/Pr1 = ',pr21

us=cs1*sqrt( ((gm-1)+(gm+1)*pr21) /(2.*gm) )
xs=x0+us*t
if not keyword_set(quiet) then print,'Vx, X of shock =',us,xs

;
; behind shock

u2=cs1*(pr21-1)*sqrt( (2/gm)/((gm+1.d0)*pr21+(gm-1.d0)) )
ro2=ro1*us/(us-u2)
pr2=pr21*pr1
if not keyword_set(quiet) then print,'Ro2,Pr2,Vx2 = ',ro2,pr2,u2


;
; contact discontinuity

uc=u2
xc=x0+uc*t
if not keyword_set(quiet) then print,'Vx, X of cont. disc. =',uc,xc

;
; behind expansion fan

u3=u2
pr3=pr2
ro3=ro4*(pr3/pr4)^(1./gm)
cs3=sqrt(gm*pr3/ro3)
if not keyword_set(quiet) then print,'Ro3,Pr3,Vx3 = ',ro3,pr3,u3

;
; expansion fan

phihead=u4-cs4
phitail=u3-cs3

xfanhead=x0+phihead*t
xfantail=x0+phitail*t
if not keyword_set(quiet) then print,'X of exp. fan =',xfanhead,xfantail

whr=where( (x gt xfanhead) and (x lt xfantail), cnt )
if (cnt le 0) then goto,end_of_expfan

  phi=(x(whr)-x0)/t
  ;xfan=x0+phi*t
  csfan=(gm-1)/(gm+1) * (2.*cs4/(gm-1) -phi)
  rofan= (csfan^2/gm)^(1./(gm-1))
  prfan= (csfan^2/gm)^(gm/(gm-1))
  ufan= csfan + phi

  ro(whr)=rofan
  pr(whr)=prfan
  vx(whr)=ufan

end_of_expfan:

;
; input values

whr=where( (x le xfanhead), cnt)
if (cnt ge 1) then begin
  ro(whr)=ro4
  pr(whr)=pr4
  vx(whr)=u4
endif

whr=where( (x ge xfantail) and (x le xc), cnt)
if (cnt ge 1) then begin
  ro(whr)=ro3
  pr(whr)=pr3
  vx(whr)=u3
endif

whr=where( (x ge xc) and (x le xs), cnt)
if (cnt ge 1) then begin
  ro(whr)=ro2
  pr(whr)=pr2
  vx(whr)=u2
endif

whr=where( (x ge xs), cnt)
if (cnt ge 1) then begin
  ro(whr)=ro1
  pr(whr)=pr1
  vx(whr)=u1
endif

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;gm=1.4d0
;ro1=0.125d0
;pr1=0.1d0
;x0=0.
;t=0.14154
;
;ix=1001
;x=1./(ix-1)*findgen(ix)
;x=x-max(x)/2.
;
;shktb_analytic,roan,pran,vxan,t,x0,x,gm,ro1,pr1,/quiet
;
;end
