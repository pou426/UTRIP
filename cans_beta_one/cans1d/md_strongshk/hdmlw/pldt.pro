
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
nt0=38<(n_elements(t)-1)

plot,x,vx(*,0),yr=[-10,25],title='vx',linest=1
oplot,x,vx(*,nt0);,psym=4

plot,x,ro(*,0),title='De',yrange=[0.,7],linest=1,/yst
oplot,x,ro(*,nt0);,psym=4

;plot,x,pr(*,0),title='Pr',/ylog,yrange=[10,10000],linest=1
;oplot,x,pr(*,nt0),psym=4

;plot,x,te(*,0),/ylog,yrange=[10.,1.e5],title='Te',linest=1
;oplot,x,te(*,nt0),psym=4

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
