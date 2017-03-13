
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
!p.multi=[0,2,2]
!x.style=1
nt0=n_elements(t)-1

plot,x,pr[*,0],title='Pr',yr=[0.5,0.7]
for n=0,nx-1 do oplot,x,pr[*,n]

plot,x,ro[*,0],title='De',yr=[0.9,1.1]
for n=0,nx-1 do oplot,x,ro[*,n]

plot,x,te[*,0],title='Te',yr=[0.9,1.1]
for n=0,nx-1 do oplot,x,te[*,n]

plot,x,vx[*,0],title='vx',yr=[-0.1,0.1]
for n=0,nx-1 do oplot,x,vx[*,n]

!p.multi=0
!x.style=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
