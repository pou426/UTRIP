
if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
x0=-1.0*sin(thini)*cos(phini) & x1= 1.0*sin(thini)*cos(phini)
y0=-1.0*sin(thini)*sin(phini) & y1= 1.0*sin(thini)*sin(phini)
z0=-1.0*cos(thini)            & z1= 1.0*cos(thini)

xi=x0+(x1-x0)/22.*findgen(23)
yi=y0+(y1-y0)/22.*findgen(23)
zi=z0+(z1-z0)/22.*findgen(23)

s=sqrt((xi-x0)^2+(yi-y0)^2+(zi-z0)^2)
s=s-max(s)/2

mx=n_elements(s)

ires=interpol(findgen(ix),x,xi)
jres=interpol(findgen(jx),y,yi)
kres=interpol(findgen(kx),z,zi)

ro_s=interpolate(ro,ires,jres,kres)
scl=fltarr(nx)
;scl(1:*)=t(1:*)^0.5  ; uniform conductivity
scl(1:*)=t(1:*)^(2./9.)  ; Spitzer conductivity

for n=0,nx-1 do begin
  pr_s=interpolate(pr(*,*,*,n),ires,jres,kres)
  te_sl0=gm*pr_s/ro_s
  if (n eq 0) then te_sl=fltarr(mx,nx)
  te_sl(*,n)=te_sl0
endfor

x_s=fltarr(mx,nx)
for n=1,nx-1 do x_s(*,n)=s/scl(n)

te_s=fltarr(mx,nx)
for n=1,nx-1 do te_s(*,n)=te_sl(*,n)*scl(n)

!p.multi=[0,1,2]
!x.style=1
nt0=1
nt1=nx

plot,s,te_sl(*,0) ,yrange=[0.,1.2],title='Te',linest=1,xrange=[0,1]
for n=nt0,nt1-1 do oplot,s,te_sl(*,n)

plot,x_s(*,1),te_s(*,1) ,title='Te-scaled',linest=1,yrange=[0,1.0],xrange=[0,1.6]
for n=nt0,nt1-1 do oplot,x_s(*,n),te_s(*,n)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  filepng='pldt1d.png'
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

!x.style=0
!p.multi=0

end
