pro slice_1dfrom3d,data1d,s,data3d,xyz0,xyz1,mx=mx,x=x,y=y,z=z

sz=size(data3d)

case sz(0) of
  3: begin ; single time step
    nx=1
  end
  4: begin ; multiple time steps
    nx=sz[4]
  end
  else: begin ; illeagal
    print,' Error: The input data array is not 2-dimensional'
    return
  end
endcase

ix=sz[1] & jx=sz[2] & kx=sz[3]

mx0=max(sz(1:3))
if (n_elements(mx) eq 0) then mx=mx0

if (n_elements(x) eq 0) then x=fltarr(ix)
if (n_elements(y) eq 0) then y=fltarr(jx)
if (n_elements(z) eq 0) then z=fltarr(kx)

x0=xyz0[0] & y0=xyz0[1] & z0=xyz0[2]
x1=xyz1[0] & y1=xyz1[1] & z1=xyz1[2]


xi=x0+(x1-x0)/float(mx-1)*findgen(mx)
yi=y0+(y1-y0)/float(mx-1)*findgen(mx)
zi=z0+(z1-z0)/float(mx-1)*findgen(mx)


s=sqrt((xi-x0)^2+(yi-y0)^2+(zi-z0)^2)

ires=interpol(findgen(ix),x,xi)
jres=interpol(findgen(jx),y,yi)
kres=interpol(findgen(kx),z,zi)

data1d=fltarr(mx,nx)

for n=0,nx-1 do begin
  data1d0=interpolate(data3d[*,*,*,n],ires,jres,kres)
  data1d(*,n)=data1d0
endfor

data1d=reform(data1d)

return
end
