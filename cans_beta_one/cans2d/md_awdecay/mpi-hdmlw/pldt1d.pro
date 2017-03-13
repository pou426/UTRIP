
if (n_elements(plot_choice) eq 0) then plot_choice='x'
if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  filepng='pldt1d.png'
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
x0=min(x)*cos(thini) & x1= min(x)*cos(thini)
y0=min(y)*sin(thini) & y1= max(y)*sin(thini)

slice_1dfrom2d,ro1d,s,ro,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,vx1d,s,vx,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,vy1d,s,vy,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,vz1d,s,vz,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,bx1d,s,bx,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,by1d,s,by,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,bz1d,s,bz,[x0,y0],[x1,y1],x=x,y=y

vs1d=vx1d*cos(thini)+vy1d*sin(thini)
vt1d=vz1d
bs1d=bx1d*cos(thini)+by1d*sin(thini)
bt1d=bz1d


!p.multi=[0,2,2]
!x.style=1

n=100

plot,s,ro1d(*,n),title='De'
oplot,s,ro1d(*,0),linest=1

plot,s,vs1d(*,0),title='Vs',linest=1,yrange=2.d-1*[-1,1]
oplot,s,vs1d(*,n)

plot,s,vt1d(*,0),title='Vt',linest=1
oplot,s,vt1d(*,n)

plot,s,bt1d(*,0)/sqrt(!pi*4),title='Bt/sqrt(4pi)',linest=1
oplot,s,bt1d(*,n)/sqrt(!pi*4)

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
