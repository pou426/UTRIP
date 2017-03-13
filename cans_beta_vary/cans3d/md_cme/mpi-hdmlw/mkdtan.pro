ro_an=fltarr(ix,jx,nx)
pr_an=fltarr(ix,jx,nx)
vx_an=fltarr(ix,jx,nx)
vy_an=fltarr(ix,jx,nx)
vz_an=fltarr(ix,jx,nx)
bx_an=fltarr(ix,jx,nx)
by_an=fltarr(ix,jx,nx)
bz_an=fltarr(ix,jx,nx)
az_an=fltarr(ix,jx,nx)

for n=0,nx-1 do begin
  time0=t(n)
  cme_analytic,x,y,time0,ro0,pr0,vx0,vy0,vz0,bx0,by0,bz0,az0 $
   ,zetac_0=zetac_0,rnu_0=rnu_0,eta_0=eta_0,a0_0=a0_0,rlambda_0=rlambda_0 $
   ,rMsun=rMsun,rr0=rr0,rn0_0=rn0_0

  ro_an(*,*,n)=ro0
  pr_an(*,*,n)=pr0
  vx_an(*,*,n)=vx0
  vy_an(*,*,n)=vy0
  vz_an(*,*,n)=vz0
  bx_an(*,*,n)=bx0
  by_an(*,*,n)=by0
  bz_an(*,*,n)=bz0
  az_an(*,*,n)=az0
endfor
end
