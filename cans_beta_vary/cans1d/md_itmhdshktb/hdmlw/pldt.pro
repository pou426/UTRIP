
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
;!p.charsize=1.5
!x.style=1

nt0=n_elements(t)-1
;nt0=5
nt1=nt0+1
n=n_elements(t)-1

plot,x,ro(*,0),title='De',yrange=[0,1.2],linest=1
;for n=nt0,nt1-1 do oplot,x,ro(*,n)
oplot,x,ro(*,n)

plot,x,vx(*,0),yr=[-1.2,1.2],title='Vx',linest=1,/yst
;for n=nt0,nt1-1 do oplot,x,vx(*,n)
oplot,x,vx(*,n)

plot,x,vy(*,0),yr=[-2.5,0.2],title='Vy',linest=1,/yst
;for n=nt0,nt1-1 do oplot,x,vy(*,n)
oplot,x,vy(*,n)

plot,x,by(*,0)/sqrt(!pi*4),yr=[-1.5,1.5],title='By/sqrt(4pi)',linest=1,/yst
;for n=nt0,nt1-1 do oplot,x,by(*,n)/sqrt(!pi*4)
oplot,x,by(*,n)/sqrt(!pi*4)

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

end
