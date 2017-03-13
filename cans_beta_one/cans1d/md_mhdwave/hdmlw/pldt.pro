
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
!p.multi=[0,3,2]
!p.charsize=2
!x.style=1
!x.margin=[5,2]
n=n_elements(t)-1

plot,x,pr(*,0),title='Pr',yrange=[0.55,0.65],linest=1
;for n=0,nx-1 do oplot,x,pr(*,n)
oplot,x,pr(*,nx-1)

plot,x,ro(*,0),title='De',yrange=[0.9,1.1],linest=1
;for n=0,nx-1 do oplot,x,ro(*,n)
oplot,x,ro(*,nx-1)

plot,x,te(*,0) ,yrange=[0.9,1.1],title='Te',linest=1
;for n=0,nx-1 do oplot,x,te(*,n)
oplot,x,te(*,nx-1)

plot,x,vx(*,0),yr=[-0.03,0.03],title='Vx',linest=1,/yst
;for n=0,nx-1 do oplot,x,vx(*,n)
oplot,x,vx(*,nx-1)

plot,x,vy(*,0),yr=[-0.03,0.03],title='Vy',linest=1,/yst
;for n=0,nx-1 do oplot,x,vy(*,n)
oplot,x,vy(*,nx-1)

plot,x,by(*,0)/sqrt(!pi*4),yr=0.78+[-0.03,0.03],title='By/sqrt(4pi)',linest=1,/yst
;for n=0,nx-1 do oplot,x,by(*,n)/sqrt(!pi*4)
oplot,x,by(*,nx-1)/sqrt(!pi*4)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
!p.charsize=1
!x.style=0
!x.margin=0

end
