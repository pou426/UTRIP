function sedov_analytic_func,x,y

common sedov_analytic_param,gmc

gm=gmc

  v=y[0]
  omega=y[1]
  pi=y[2]

b=[ -v*(v-1.d0)-2.d0*pi/omega, -3.d0*v, -2.d0*(v-1.d0)]

a=[[v-2.d0/5.d0, 0.d0, 1.d0/omega], $
   [1.d0, (v-2.d0/5.d0)/omega, 0.d0], $
   [0.d0, -(v-2.d0/5.d0)*gm/omega, (v-2.d0/5.d0)/pi]]

dydx=reform(invert(a)##b)/x

return,dydx

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sedov_analytic,omega,pi,v,xi,gm

common sedov_analytic_param,gmc

gmc=gm

; Normalized solutions

x0=1.d0
y0=[2.d0/(gm+1.d0)*2.d0/5.d0, (gm+1.d0)/(gm-1.d0), 2.d0/(gm+1)*(2.d0/5.d0)^2]

nx=1000L
h=-0.001d0

y0ar=fltarr(3,nx)
x0ar=fltarr(nx)

x0ar[0]=x0
y0ar[*,0]=y0

for n=1L,nx-1 do begin
  dydx=sedov_analytic_func(x0,y0)
  result=rk4(y0,dydx,x0,h,'sedov_analytic_func',/double)
  x0=x0+h
  y0=result
  x0ar[n]=x0
  y0ar[*,n]=y0
endfor

lambda=reverse(x0ar)
v=reverse(reform(y0ar[0,*]))
omega=reverse(reform(y0ar[1,*]))
pi=reverse(reform(y0ar[2,*]))

;z=omega^(1-gm)*pi*(v-0.4)*lambda^5

xi0=(1/(total((0.5*omega*v^2+pi/(gm-1.))*4.*!pi*lambda^4)*abs(h)))^0.2

xi=lambda*xi0


end
