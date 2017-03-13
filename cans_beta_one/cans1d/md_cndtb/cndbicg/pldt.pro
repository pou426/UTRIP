
if (n_elements(plot_choice) eq 0) then plot_choice='x'
if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  filepng='pldt.png'
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
scl=fltarr(nx)
;scl[1:*]=t[1:*]^0.5  ; uniform conductivity
scl[1:*]=t[1:*]^(2./9.)  ; Spitzer conductivity

x_s=fltarr(ix,nx)
for n=1,nx-1 do x_s[*,n]=x/scl[n]

te_s=fltarr(ix,nx)
for n=1,nx-1 do te_s[*,n]=te[*,n]*scl[n]


!x.style=1
nt0=1
nt1=nx

!p.multi=[0,1,2]

plot,x,te[*,0] ,yrange=[0.,1.2],title='Te',linest=1
for n=nt0,nt1-1 do oplot,x,te[*,n]

plot,x_s[*,1],te_s[*,1] ,title='Te-scaled',linest=1,yrange=[0,1.0]
for n=nt0,nt1-1 do oplot,x_s[*,n],te_s[*,n]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
end
