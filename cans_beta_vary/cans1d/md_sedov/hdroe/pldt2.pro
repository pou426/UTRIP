
if (n_elements(plot_choice) eq 0) then plot_choice='x'
if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  filepng='pldt2.png'
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

n=6

ro0=1.

enttl=(sqrt(!pi)*wexp)^3/(gm-1)

sedov_analytic,omega,pi,v,xi,gm

xan=xi*t[n]^(0.4)*(enttl/ro0)^(0.2)
roan=ro0*omega
pran=ro0*pi*xan^2/t[n]^2
vxan=v*xan/t[n]

!p.multi=[0,1,3]

!p.charsize=2
!x.range=[0,1]

plot,x,ro[*,n]
oplot,xan,roan,linest=2
plot,x,pr[*,n]
oplot,xan,pran,linest=2
plot,x,vx[*,n]
oplot,xan,vxan,linest=2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
!p.charsize=1
!x.range=0

end
