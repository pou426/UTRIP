
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

nt=nx-1

plot,x,ro(*,0),/xlog,/ylog,xr=[10,1000],yr=10.^[-9,0],title='ro',linest=1
for n=1,nx-1 do oplot,x,ro(*,n)
oplot,x,ro(*,nt)
oplot,x,0.03*x^(-2),linest=2

plot,x,vx(*,0),yr=[0,2],title='Vx',linest=1,/xst
for n=1,nx-1 do oplot,x,vx(*,n)
oplot,x,vx(*,nt)
oplot,x,sqrt(te(*,nt)),linest=2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
end
