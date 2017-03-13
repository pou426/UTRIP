
if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (plot_choice eq 'x') then plot_interactive=1 else plot_interactive=0

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
j=4

enttl=(sqrt(!pi)*wexp)^3/(gm-1)
scl=fltarr(nx)
scl(1:*)=1.12*enttl^0.2*t(1:*)^0.4

d_s=fltarr(nx)
d_s(1:*)=0.4*scl(1:*)/t(1:*)

x_s=fltarr(ix,nx)
for n=1,nx-1 do x_s(*,n)=x/scl(n)

pr_s=pr
ro_s=ro
vx_s=vx

for n=1,nx-1 do pr_s(*,n)=pr(*,j,n)/2*(gm+1)/d_s(n)^2
for n=1,nx-1 do ro_s(*,n)=ro(*,j,n)/(gm+1)*(gm-1)
for n=1,nx-1 do vx_s(*,n)=vx(*,j,n)/2*(gm+1)/d_s(n)

!p.multi=[0,3,2]
;!x.range=[0.,0.2]
!x.style=1
!p.charsize=1.5

nt0=1
nt1=n_elements(t)
dn=1

!x.range=[0,0.4]
plot,x,pr(*,j,0),title='Pr',yrange=[1.e-4,1],linest=1,/yl
for n=nt0,nt1-1,dn do oplot,x,pr(*,j,n)

plot,x,ro(*,j,0),title='De',yrange=[0,5],linest=1
for n=nt0,nt1-1,dn do oplot,x,ro(*,j,n)

plot,x,vx(*,j,0),yr=[-0.05,0.15],title='vx',linest=1
for n=nt0,nt1-1,dn do oplot,x,vx(*,j,n)

!x.range=[0,1.5]
plot,x_s(*,1),pr_s(*,1),title='Pr-scaled',yrange=[0.,2.0],linest=1
for n=nt0,nt1-1,dn do oplot,x_s(*,n),pr_s(*,n)

plot,x_s(*,1),ro_s(*,1),title='De-scaled',yrange=[0,1.2],linest=1
for n=nt0,nt1-1,dn do oplot,x_s(*,n),ro_s(*,n)

plot,x_s(*,1),vx_s(*,1),title='Vx-scaled',yr=[0.,2.0],linest=1
for n=nt0,nt1-1,dn do oplot,x_s(*,n),vx_s(*,n)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  filepng='idl.png'
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
!x.style=0
!x.range=0
!p.charsize=1.

end
