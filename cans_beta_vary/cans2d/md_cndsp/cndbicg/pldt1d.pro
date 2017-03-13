if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xy0=[0,0]
xy1=[1,1]

slice_1dfrom2d,te1d,s,te,xy0,xy1,x=x,y=z
mx=n_elements(s)

scl=fltarr(nx)
;scl(1:*)=t(1:*)^0.5  ; uniform conductivity
scl(1:*)=t(1:*)^(2./19.)  ; Spitzer conductivity

x_s=fltarr(mx,nx)
for n=1,nx-1 do x_s(*,n)=s/scl(n)

te_s=fltarr(mx,nx)
for n=1,nx-1 do te_s(*,n)=te1d(*,n)*scl(n)^3


!x.style=1

!p.multi=[0,1,2]

plot,s,te1d(*,0) ,yrange=[0.,1],title='Te',linest=1
for n=1,nx-1 do oplot,s,te1d(*,n)

plot,x_s(*,1),te_s(*,1) ,title='Te-scaled',linest=1,yrange=[0,0.4]
for n=1,nx-1 do oplot,x_s(*,n),te_s(*,n)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  filepng='pldt1d.png'
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
end
