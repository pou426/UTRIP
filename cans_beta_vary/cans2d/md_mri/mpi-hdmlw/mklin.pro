kzar=[1,2,3,4,5,6]/(z(jx-1)-z(0))
mx=n_elements(kzar)
ampar=fltarr(nx,mx)

dz=z-shift(z,1) & dz(0)=dz(1)
dx=x-shift(x,1) & dx(0)=dx(1)

for m=0,mx-1 do begin

  kz=kzar(m)
  amp=fltarr(nx)

  for n=0,nx-1 do begin

    f=fltarr(ix)

    for i=0,ix-1 do begin
      fcos=total(bx(i,*,n)*cos(2*!pi*kz*z)*dz)
      fsin=total(bx(i,*,n)*sin(2*!pi*kz*z)*dz)
      f(i)=sqrt(fcos^2+fsin^2)
    endfor

    amp(n)=total(f*dx)/(x(ix-1)-x(0))

  endfor

  ampar(*,m)=amp

endfor



end
