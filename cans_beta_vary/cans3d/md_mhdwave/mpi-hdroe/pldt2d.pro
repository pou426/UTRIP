
if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
unitx=sin(thini)*cos(phini)
unity=sin(thini)*sin(phini)
unitz=cos(thini)

ebx0=1.
eby0=2.
ebz0=-(ebx0*unitx+eby0*unity)/unitz
ebb0=sqrt(ebx0^2+eby0^2+ebz0^2)
ebx=ebx0/ebb0
eby=eby0/ebb0
ebz=ebz0/ebb0

nlevels=32
color_index=(256/nlevels)*indgen(nlevels)

dmin=0.9995d0 & dmax=1.0005d0
dlevel=(dmax-dmin)/nlevels
levels1=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)

xrot=-asin(eby)
yrot= acos(eby)
zrot=0.

slice=fltarr(ix,jx,nx)

for n=0,nx-1 do begin
  data=reform(ro[*,*,*,n])
  slice(*,*,n)=extract_slice(data,ix,jx,ix/2,jx/2,kx/2,xrot,yrot,zrot)
endfor

!p.multi=[0,3,2]

;loadct,5,/silent

for m=0,5 do begin
  contour,slice[*,*,m],/fill,/xst,/yst,levels=levels1,c_colors=color_index
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  filepng='pldt2d.png'
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

end
