
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
kzero=x2i(z,0.)

!p.multi=[0,2,2]

plot,x,pr(*,jzero,kzero,0) ,yrange=[0.1,100],title='Te',linest=1,/yl,xtitle='x'
for n=nt0,nt1-1 do oplot,x,pr(*,jzero,kzero,n)

plot,y,pr(izero,*,kzero,0) ,yrange=[0.1,100],title='Te',linest=1,/yl,xtitle='y'
for n=nt0,nt1-1 do oplot,y,pr(izero,*,kzero,n)

plot,z,pr(izero,jzero,*,0) ,yrange=[0.1,100],title='Te',linest=1,/yl,xtitle='z'
for n=nt0,nt1-1 do oplot,z,pr(izero,jzero,*,n)

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
