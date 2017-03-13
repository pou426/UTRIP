kzar=[1,2,3,4,5,6]/zrgn
mx=n_elements(kzar)
ampar=fltarr(nx,mx)

dz=z-shift(z,1) & dz[0]=dz[1]
dy=y-shift(y,1) & dy[0]=dy[1]
dx=x-shift(x,1) & dx[0]=dx[1]

for m=0,mx-1 do begin

  kz=kzar(m)
  amp=fltarr(nx)

  for n=0,nx-1 do begin

    f=fltarr(ix,jx)

    for j=0,jx-1 do begin
    for i=0,ix-1 do begin
      fcos=total(bx[i,j,*,n]*cos(2*!pi*kz*z)*dz)
      fsin=total(bx[i,j,*,n]*sin(2*!pi*kz*z)*dz)
      f[i,j]=sqrt(fcos^2+fsin^2)*dx[i]*dy[j]/(xrgn*yrgn)
    endfor
    endfor

    amp(n)=total(f)

  endfor

  ampar(*,m)=amp

endfor



end
