
if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!x.style=1
nt0=1
nt1=nx

izero=x2i(x,0.)
jzero=x2i(y,0.)

!p.multi=[0,2,1]

plot,x,te(*,jzero,0) ,yrange=[0.1,100],title='Te',linest=1,/yl
for n=nt0,nt1-1 do oplot,x,te(*,jzero,n)

plot,y,te(izero,*,0) ,yrange=[0.1,100],title='Te',linest=1,/yl
for n=nt0,nt1-1 do oplot,y,te(izero,*,n)

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
