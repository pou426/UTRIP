phi=fltarr(nx)
for n=0,nx-1 do begin
  phi(n)=total(abs(by(3:ix-3,4,n)))
endfor

end
