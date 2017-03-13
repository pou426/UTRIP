
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

plot,x,pr(*,0),linest=1,yr=[0.,1.5],title='Pr'
oplot,x,pr(*,nt0/4)
oplot,x,pr(*,nt0/2)
oplot,x,pr(*,nt0)

plot,x,ro(*,0),linest=1,yr=[0.,1.5],title='De'
oplot,x,ro(*,nt0/4)
oplot,x,ro(*,nt0/2)
oplot,x,ro(*,nt0)

plot,x,te(*,0),linest=1,yr=[0,2],title='Te'
oplot,x,te(*,nt0/4)
oplot,x,te(*,nt0/2)
oplot,x,te(*,nt0)

plot,x,vx(*,0),linest=1,title='vx'
oplot,x,vx(*,nt0/4)
oplot,x,vx(*,nt0/2)
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
