mass=fltarr(nx)
dx=mkdelta(x)
dz=mkdelta(z)
dv=(x*dx)#dz
for n=0,nx-1 do begin
; mass(n)=total(reform(ro(*,*,n))*dv)
  en=reform(pr(*,*,n)/(gm-1.d0) $
     +0.5d0*ro(*,*,n)*(vx(*,*,n)^2+vz(*,*,n)^2))
  mass(n)=total(en(2:*,2:*)*dv(2:*,2:*))
endfor
plot,t,(mass-mass(0))/mass(0)
end
