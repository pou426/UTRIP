
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
xy0=[0,0]
xy1=[1,1]

slice_1dfrom2d,ro_i,s,ro,xy0,xy1,x=x,y=z
slice_1dfrom2d,pr_i,s,pr,xy0,xy1,x=x,y=z
slice_1dfrom2d,vx_i,s,vx,xy0,xy1,x=x,y=z
slice_1dfrom2d,vz_i,s,vz,xy0,xy1,x=x,y=z
vs_i=vx_i*cos(!pi/4)+vz_i*sin(!pi/4)
mx=n_elements(s)

enttl=(sqrt(!pi)*wexp)^3/(gm-1)
scl=fltarr(nx)
scl(1:*)=1.12*enttl^0.2*t(1:*)^0.4

d_s=fltarr(nx)
d_s(1:*)=0.4*scl(1:*)/t(1:*)

x_s=fltarr(mx,nx)
for n=1,nx-1 do x_s(*,n)=s/scl(n)

pr_s=pr_i
ro_s=ro_i
vs_s=vs_i

for n=1,nx-1 do pr_s(*,n)=pr_i(*,n)/2*(gm+1)/d_s(n)^2
for n=1,nx-1 do ro_s(*,n)=ro_i(*,n)/(gm+1)*(gm-1)
for n=1,nx-1 do vs_s(*,n)=vs_i(*,n)/2*(gm+1)/d_s(n)

!p.multi=[0,3,2]
;!x.range=[0.,0.2]
!x.style=1
nt0=1
nt1=n_elements(t)
dn=1

!x.range=[0,0.4]
plot,s,pr_i(*,0),title='Pr',yrange=[1.e-4,1],linest=1,/yl
for n=nt0,nt1-1,dn do oplot,s,pr_i(*,n)

plot,s,ro_i(*,0),title='De',yrange=[0,5],linest=1
for n=nt0,nt1-1,dn do oplot,s,ro_i(*,n)

plot,s,vs_i(*,0),yr=[-0.05,0.1],title='vx',linest=1
for n=nt0,nt1-1,dn do oplot,s,vs_i(*,n)

!x.range=[0,1.5]
plot,x_s(*,1),pr_s(*,1),title='Pr-scaled',yrange=[0.,1.2],linest=1
for n=nt0,nt1-1,dn do oplot,x_s(*,n),pr_s(*,n)

plot,x_s(*,1),ro_s(*,1),title='De-scaled',yrange=[0,1.2],linest=1
for n=nt0,nt1-1,dn do oplot,x_s(*,n),ro_s(*,n)

plot,x_s(*,1),vs_s(*,1),title='Vx-scaled',yr=[0.,1.2],linest=1
for n=nt0,nt1-1,dn do oplot,x_s(*,n),vs_s(*,n)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
!x.style=0
!x.range=0

end
