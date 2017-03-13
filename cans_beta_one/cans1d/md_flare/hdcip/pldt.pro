
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
nt0=0
nt1=nx
nts=1

if 1 then begin
plot,x,pr(*,0),/ylog,title='Pr'
for n=nt0,nt1-1 do oplot,x,pr(*,n)

plot,x,ro(*,0),/yl,title='De'
for n=nt0,nt1-1 do oplot,x,ro(*,n)

plot,x,te(*,0) ,/yl,yrange=[0.1,10000],title='Te'
for n=nt0,nt1-1 do oplot,x,te(*,n)

plot,x,vx(*,0),yr=[-10,20],title='vx'
for n=nt0,nt1-1 do oplot,x,vx(*,n)
endif

if 0 then begin
nt00=20
!x.range=[1000,10000]
prnml=ronml*tenml*1.e-16
venml=sqrt(tenml*1.e8)/1.e5

plot,x*rlnml*1.e-5,te(*,nt00)*tenml,title='Te',/xlog,/ylog,yr=[1.e3,1.e7]
for n=nt0,nt1-1,nts do oplot,x*rlnml*1.e-5,te(*,n)*tenml

plot,x*rlnml*1.e-5,ro(*,nt00)*ronml,/ylog,title='De',/xlog
for n=nt0,nt1-1,nts do oplot,x*rlnml*1.e-5,ro(*,n)*ronml

plot,x*rlnml*1.e-5,vx(*,nt00)*venml,yr=[-300,300],title='vx'
for n=nt0,nt1-1,nts do oplot,x*rlnml*1.e-5,vx(*,n)*venml

plot,x*rlnml*1.e-5,pr(*,nt00)*prnml,/ylog,title='Pr',yr=[10,1000]
for n=nt0,nt1-1,nts do oplot,x*rlnml*1.e-5,pr(*,n)*prnml
!x.range=0
endif

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
