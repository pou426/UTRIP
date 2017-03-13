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
!p.multi=[0,3,1]
!x.style=1
nt0=n_elements(t)
!x.range=[1.e-6,1.e2]
!p.charsize=2

plot,x(2:*),ro(2:*,0),title='De',yrange=[1.e-2,1.e7] $
 ,/xlog,/ylog,/yst
for n=0,nt0-1,1 do oplot,x(2:*),ro(2:*,n)
oplot,x(2:*),ro(2:*,x2i(t,0.293))
oplot,x(2:*),ro(2:*,x2i(t,0.326))
oplot,x(2:*),ro(2:*,x2i(t,0.327))


plot,x(2:*),-vx(2:*,0)>(1.e-4),yr=[1.e-4,1.e4],title='-vx',linest=1,/xlog,/ylog
for n=0,nt0-1,1 do oplot,x(2:*),-vx(2:*,n)>(1.e-4)
oplot,x(2:*),-vx(2:*,x2i(t,0.293))
oplot,x(2:*),-vx(2:*,x2i(t,0.326))
oplot,x(2:*),-vx(2:*,x2i(t,0.327))

plot,x(2:*),-gx(2:*,1),/xlog,/ylog,yr=[1.e-4,1.e8] $
  ,title='-gx'
for n=0,nt0-1,1 do oplot,x(2:*),-gx(2:*,n)
oplot,x(2:*),-gx(2:*,x2i(t,0.293))
oplot,x(2:*),-gx(2:*,x2i(t,0.326))
oplot,x(2:*),-gx(2:*,x2i(t,0.327))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
!x.range=0
!x.style=0
!p.charsize=1

end
