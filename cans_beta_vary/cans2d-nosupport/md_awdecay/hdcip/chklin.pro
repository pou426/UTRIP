if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

!p.charsize=2

vs=vx*cos(thini)+vy*sin(thini)

vxmax=fltarr(nx)
for n=0,nx-1 do vxmax(n)=sqrt(total(vs(*,*,n)^2))/(ix*jx)
plot,t(1:nx-2),vxmax(1:nx-2)>(1.e-16),/yl,yr=10.^[-10,0],xr=[0,3]
omega_growth=0.30 ; for k=2, B^2/B_0^2=0.9 & (Cs/CA)^2=0.1 (Goldstein 1978)
ca=sqrt(cs2/csca2)
omega_a=2*!pi*ca

oplot,t,1.e-4*exp(omega_growth*omega_a*t),linest=2
oplot,t,1.e-6*exp(omega_growth*omega_a*t),linest=2
oplot,t,1.e-8*exp(omega_growth*omega_a*t),linest=2

!p.charsize=1
!x.style=0

if (plot_choice eq 'png') then begin
  filepng='chklin.png'
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
end
