
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

slice_1dfrom3d,te1d,s,te,[x0,y0,z0],[x1,y1,z1],x=x,y=y,z=z
mx=n_elements(s)
s=s-max(s)/2.

scl=fltarr(nx)
;scl(1:*)=t(1:*)^0.5  ; uniform conductivity
scl(1:*)=t(1:*)^(2./9.)  ; Spitzer conductivity

x_s=fltarr(mx,nx)
for n=1,nx-1 do x_s(*,n)=s/scl(n)

te_s=fltarr(mx,nx)
for n=1,nx-1 do te_s(*,n)=te1d(*,n)*scl(n)

!p.multi=[0,1,2]
!x.style=1
nt0=1
nt1=nx

plot,s,te1d(*,0) ,yrange=[0.,1.2],title='Te',linest=1,xrange=[0,1]
for n=nt0,nt1-1 do oplot,s,te1d(*,n)

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
