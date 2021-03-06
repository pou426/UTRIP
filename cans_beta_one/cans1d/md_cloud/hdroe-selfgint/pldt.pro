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
!p.multi=[0,2,1]
!x.style=1
nt0=n_elements(t)

plot,x[ix/2:*],ro[ix/2:*,0],title='De' $
  ,yrange=[0.1,1.e2],linest=1,/xlog,/ylog
for n=0,nt0-1,1 do oplot,x[ix/2:*],ro[ix/2:*,n]


plot,x[ix/2:*],-vx[ix/2:*,0],title='-vx' $
,linest=1,/xlog,yrange=[-5,15]
for n=0,nt0-1,1 do oplot,x[ix/2:*],-vx[ix/2:*,n]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
!x.style=0

end
