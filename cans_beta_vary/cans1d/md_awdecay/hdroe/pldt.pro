
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
nt0=0
dn=5
nt1=n_elements(t)-1

plot,x,ro(*,0),title='De',yrange=[0,2],linest=1,/xst
for n=nt0,nt1-1,dn do oplot,x,ro(*,n)

plot,x,vx(*,0),yr=[-0.5,0.5],title='Vx',linest=1,/yst,/xst
for n=nt0,nt1-1,dn do oplot,x,vx(*,n)

plot,x,vy(*,0),yr=[-5,5],title='Vy',linest=1,/yst,/xst
for n=nt0,nt1-1,dn do oplot,x,vy(*,n)

plot,x,by(*,0)/sqrt(!pi*4),yr=[-5,5],title='By/sqrt(4pi)',linest=1,/yst,/xst
for n=nt0,nt1-1,dn do oplot,x,by(*,n)/sqrt(!pi*4)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

end
