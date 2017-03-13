
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
!P.MULTI=[0,2,2]

!x.style=1
!y.style=1
!y.range=[0,t(nx-1)*2.0]
!x.range=[0.01,1.d2]

plot, zz, vx(*,0)*50+t(0),/xlog,title='Vx'
for n=0,nx-1,1 do oplot, zz, vx(*,n)*50+t(n)

plot, zz, (alog10(ro(*,0))+7.0)*8+t(0),/xlog,title='log10(ro)'
for n=0,nx-1,1 do oplot, zz, (alog10(ro(*,n))+7.0)*8+t(n)

plot, zz, vy(*,0)*50+t(0),/xlog,title='Vy'
for n=0,nx-1,1 do oplot, zz, vy(*,n)*50+t(n)

plot, zz, -by(*,0)/bx*10+t(0),/xlog,title='-By'
for n=0,nx-1,1 do oplot, zz, -by(*,n)/bx*10+t(n)


!p.multi=0
!x.style=0
!y.style=0
!y.range=0
!x.range=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end
