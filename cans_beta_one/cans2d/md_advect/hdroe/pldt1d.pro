if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (plot_choice eq 'png') then begin
  filepng='pldt1d.png'
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

xy0=[-0.5,-0.5]
xy1=[ 0.5, 0.5]

slice_1dfrom2d,ro1d,s,ro,xy0,xy1,x=x,y=y

plot,s,ro1d(*,0),title='De',yrange=[0,1.2],linest=1
for n=0,nx-1 do oplot,s,ro1d(*,n)

if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif

end
