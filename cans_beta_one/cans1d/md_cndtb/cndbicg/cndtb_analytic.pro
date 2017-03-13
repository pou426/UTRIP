function cndtb_analytic_func,x,y

common cndtb_analytic_param,idxcndc

idxcnd=idxcndc

  f=y[0]
  g=y[1]

dydx=[g/f^idxcnd,-1/(idxcnd+2.)*(x*g/f^idxcnd+f)]

return,dydx

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro cndtb_analytic,f,xi,idxcnd

common cndtb_analytic_param,idxcndc

idxcndc=idxcnd

; Normalized solutions

y00=1.d0
dy00=-0.001

for m=0,300 do begin

x0=0.d0
y0=[y00,0.d0]

nx=999L
h=0.01d0

y0ar=fltarr(2,nx)
x0ar=fltarr(nx)

x0ar[0]=x0
y0ar[*,0]=y0

for n=1L,nx-1 do begin
  dydx=cndtb_analytic_func(x0,y0)
  result=rk4(y0,dydx,x0,h,'cndtb_analytic_func',/double)
  x0=x0+h
  y0=result
  if (not finite(y0[0])) then goto, L100
  x0ar[n]=x0
  y0ar[*,n]=y0
endfor
L100:
xi=x0ar[0:n-1]
f=reform(y0ar[0,0:n-1])

dxi=xi-shift(xi,1) & dxi[0]=dxi[1]
fttl=total(f*dxi)
;print,fttl
if (fttl le 1.d0) then goto,L200

y00=y00+dy00

endfor
L200:

end
