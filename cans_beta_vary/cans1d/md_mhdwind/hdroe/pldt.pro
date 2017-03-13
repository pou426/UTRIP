
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
!p.charsize=1.5
!x.margin=[5,3]
!y.margin=[2,2]

dn=10
nt=nx-1
;nt=50

plot, x, ro(*,nt),/ylog,title='ro',/xlog,/xst
oplot,x,x^(-2),linest=2
oplot,[6,6],10.^!y.crange,linest=1

plot, x, pr(*,nt),/ylog,title='pr',/xlog,/xst
oplot,[6,6],10.^!y.crange,linest=1
plot, x, te(*,nt),/ylog,title='te',/xlog,/xst
oplot,[6,6],10.^!y.crange,linest=1
plot, x, vx(*,nt),title='vx',/xlog,/xst
oplot,[6,6],!y.crange,linest=1
plot, x, -by(*,nt),/ylog,title='-by',/xlog,/xst
oplot,[6,6],10.^!y.crange,linest=1
oplot,x,0.3*x^(-1),linest=2
plot, x, vy(*,nt),/ylog,title='vy',/xlog,/xst
oplot,[6,6],10.^!y.crange,linest=1
oplot,x,3.*x^(-1),linest=2


if 0 then begin
plot, x, ro(*,nt),/ylog,title='ro',yr=10.^[-9,0]
for n=0,nx-1,dn do oplot, x, ro(*,n)

plot, x, pr(*,nt),/ylog,title='pr',yr=10.^[-10,0]
for n=0,nx-1,dn do oplot, x, pr(*,n)

plot, x, te(*,nt),/ylog,title='te'
for n=0,nx-1,dn do oplot, x, te(*,n)

plot, x, vx(*,nt),title='vx',yr=[-1,5]
for n=0,nx-1,dn do oplot, x, vx(*,n)

plot, x, by(*,nx-1),title='by',yr=[-3,0.]
for n=0,nx-1,dn do oplot, x, by(*,n)

plot, x, vy(*,nx-1),title='vy',yr=[-4,2]
for n=0,nx-1,dn do oplot, x, vy(*,n)
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

!p.multi=0

!p.charsize=1.0
!x.margin=[10.,3.]
!y.margin=[4.,2.]


end
