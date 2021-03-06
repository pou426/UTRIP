
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
!P.MULTI=[0,3,1]

!x.style=1
!y.style=1
!y.range=[0,t(nx-1)*1.2]
!x.range=[0,300]
!p.charsize=1.5


plot, x, (alog10(ro(*,0))+7.0)*3+t(0),title='log10(ro)'
for n=0,nx-1 do oplot, x, (alog10(ro(*,n))+7.0)*3+t(n)

plot, x, vx(*,0)+t(0),title='Vx'
for n=0,nx-1 do oplot, x, vx(*,n)+t(n)

plot, x, vy(*,0)+t(0),title='Vy'
for n=0,nx-1 do oplot, x, vy(*,n)+t(n)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

!p.multi=0
!x.style=0
!y.style=0
!y.range=0
!x.range=0


end
