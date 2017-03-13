
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

plot,x,pr(*,0),title='Pr',yrange=[0,1.2],linest=1
oplot,x,pr(*,nt0)

plot,x,ro(*,0),title='De',yrange=[0,1.2],linest=1
oplot,x,ro(*,nt0)


plot,x,te(*,0),yrange=[0.,2],title='Te',linest=1
oplot,x,te(*,nt0)

plot,x,vx(*,0),yr=[-1.0,1.5],title='vx',linest=1
oplot,x,vx(*,nt0)

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
