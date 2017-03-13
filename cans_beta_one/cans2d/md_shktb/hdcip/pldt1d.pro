
if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
shktb_analytic,tan,xan,vxan,roan,pran
tean=pran/roan*gm

x0=min(x)*cos(thini) & x1= max(x)*cos(thini)
y0=min(y)*sin(thini) & y1= max(y)*sin(thini)

slice_1dfrom2d,ro1d,s,ro,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,pr1d,s,pr,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,te1d,s,te,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,vx1d,s,vx,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,vy1d,s,vy,[x0,y0],[x1,y1],x=x,y=y
vs1d= vx1d*cos(thini)+vy1d*sin(thini)
vt1d=-vx1d*sin(thini)+vy1d*cos(thini)
s=s-max(s)/2.


!p.multi=[0,2,2]
!x.style=1

n=nx-1

plot,xan,pran,title='Pr',yrange=[0,1.2],linest=1
oplot,s,pr1d(*,n),psym=4

plot,xan,roan,title='De',yrange=[0,1.2],linest=1
oplot,s,ro1d(*,n),psym=4

plot,xan,tean ,yrange=[0.,2],title='Te',linest=1
oplot,s,te1d(*,n),psym=4

plot,xan,vxan,yr=[-0.5,1],title='vx',linest=1
oplot,s,vs1d(*,n),psym=4

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  filepng='pldt1d.png'
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
!x.style=0

end
