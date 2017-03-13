if (n_elements(plot_choice) eq 0) then plot_choice='x'
if (n_elements(plot_interactive) eq 0) then plot_interactive=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

iwx=3 & jwx=2
if (plot_interactive eq 1) then $
read,' Plot columns & rows ? : ',iwx,jwx

nstart=0
if (plot_interactive eq 1) then $
read,'  start step ? : ',nstart

nskip=3
if (plot_interactive eq 1) then $
read,'  step interval ? : ',nskip

nlevels=32
color_index=(256/nlevels)*indgen(nlevels)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pcharsize=!p.charsize
!p.charsize=2.0

; plot

ax=30 & az=30
position=[0.01,0.01,0.99,0.99,0.26,0.74]
position=[0.01,0.01,0.99,0.99,0.01,0.99]

xrange=[min(x),max(x)] & yrange=[min(y),max(y)] & zrange=[min(z),max(z)]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


if (plot_choice eq 'png') then begin
  xsize=640 & ysize=480
endif

if (plot_choice eq 'x') then begin
  dname=!d.name
  if (!d.window eq -1) then window  
  xsize=!d.x_size & ysize=!d.y_size
endif

img=bytarr(xsize,ysize)
xdw=xsize/iwx & ydw=ysize/jwx
mwx=iwx*jwx

set_plot,'z',/copy
device,set_resolution=[xdw,ydw]
device,set_character_size=[6,9]
tvlct,red,green,blue,/get

n=nstart
mw=0

while((mw lt mwx) and (n lt nx)) do begin
; erase

  ; setup 3D coordinate

  scale3,ax=ax,az=az
  surface,/t3d,/nodata,/xst,/yst,/zst,/save $
    ,[[0,0],[0,0]],[0,1],[0,1] $
    ,position=position,xrange=xrange,yrange=yrange,zrange=zrange $
    ,xtitle='X' ,ytitle='Y' ,ztitle='Z'
  plots,/t3d,!x.crange,[1,1]*yrange[1],[1,1]*zrange[0]
  plots,/t3d,!x.crange,[1,1]*yrange[1],[1,1]*zrange[1]
  plots,/t3d,[1,1]*xrange[1],!y.crange,[1,1]*zrange[0]
  plots,/t3d,[1,1]*xrange[1],!y.crange,[1,1]*zrange[1]
  plots,/t3d,[1,1]*xrange[1],[1,1]*yrange[0],!z.crange
  plots,/t3d,[1,1]*xrange[1],[1,1]*yrange[1],!z.crange


  ; B-lines
  mx=6 & radius=3
  s0=2.*!pi/mx*findgen(mx)
  y0=0. &  z0=-10.
  zplot=z0+radius*cos(s0)
  yplot=radius*sin(s0)

  xplot=[fltarr(mx)-30.,fltarr(mx)+30.]
  yplot=[yplot,yplot]
  zplot=[zplot,zplot]

  intline,bx(*,*,*,n),by(*,*,*,n),bz(*,*,*,n) $
    ,x,y,z,xplot,yplot,zplot $
    ,index,xline,yline,zline $
    ,n_points=3000
  mx=n_elements(index)
  for m=0,mx-1 do begin
    m1=index(m)
    if (m eq 0) then m0=0 else m0=index(m-1)+1
;   plots,/t3d,xline[m0:m1],yline[m0:m1],zline[m0:m1],/data
  if (m1 gt m0) then $
    plot_tube,/t3d,xline[m0:m1],yline[m0:m1],zline[m0:m1],radius=0.5
  endfor

  ; Y slice
  t3d,/yzexch
  yslice=0.
  zvalue=yslice*!y.s[1]+!y.s[0]
  zvalue=position[3]+0.01
  levels1=-5+0.35*findgen(16)
  contour,/noerase $
   ,/t3d,position=position[[0,4,2,5]],zvalue=zvalue $
   ,reform(alog10(ro(*,x2i(y,yslice),*,n)>1.e-10)),x,z  $
   ,levels=levels1,xrange=xrange,yrange=zrange $
   ,xst=5,yst=5,/cell_fill,c_colors=color_index
   fccnve,/noerase,/noaxis,/nolabels $
     ,/t3d,position=position[[0,4,2,5]],zvalue=zvalue-0.01 $
     ,nofilled=1,novector=0,nocontour=0 $
     ,xcord=x,ycord=z $
     ,xrange=xrange,yrange=zrange $
     ,reform(alog10(pm(*,x2i(y,yslice),*,n)>1.e-10)) $
     ,levels2=-5+0.5*findgen(16),c_thick=1 $
     ,reform(vx(*,x2i(y,yslice),*,n)),reform(vz(*,x2i(y,yslice),*,n)) $
     ,iskip=1,jskip=5,scale=5,limit=0.5,v_thick=1.5,v_color=0
  t3d,/yzexch

  ; X slice
  t3d,/xzexch
  xslice=0
  zvalue=xslice*!x.s[1]+!x.s[0]
  zvalue=position[2]+0.01
  levels1=0.01*findgen(16)
  contour,/noerase $
   ,/t3d,position=position[[4,1,5,3]] $
   ,transpose(bx(x2i(x,xslice),*,*,n)>0)>levels1[1],z,y  $
   ,levels=levels1,xrange=zrange,yrange=yrange $
   ,xst=5,yst=5,/cell_fill,c_colors=color_index,zvalue=zvalue
   fccnve,/noerase,/noaxis,/nolabels $
     ,/t3d,position=position[[4,1,5,3]],zvalue=zvalue-0.01 $
     ,nofilled=1,novector=0,nocontour=1 $
     ,xcord=z,ycord=y $
     ,xrange=zrange,yrange=yrange $
     ,transpose(reform(vz(x2i(x,xslice),*,*,n))) $
     ,transpose(reform(vy(x2i(x,xslice),*,*,n))) $
     ,iskip=5,jskip=4,scale=5,limit=0.5,v_thick=1.5,v_color=255
  t3d,/xzexch

  ; Z slice
  zslice=0
  zvalue=zslice*!z.s[1]+!z.s[0]
  zvalue=position[4]-0.01
  levels1=-0.41+0.05*findgen(15)
  contour,/noerase $
   ,/t3d,position=position[[0,1,2,3]] $
   ,bz(*,*,x2i(z,zslice),n)>levels1[1],x,y  $
   ,levels=levels1,xrange=xrange,yrange=yrange $
   ,xst=5,yst=5,/cell_fill,c_colors=color_index,zvalue=zvalue

  put_time,t[n],charsize=1.5,xp=0.1,yp=0.9

  d=tvrd()

  jw=mw/iwx
  iw=mw-jw*iwx
  img[iw*xdw:(iw+1)*xdw-1,(jwx-1-jw)*ydw:(jwx-jw)*ydw-1]=d

   n=n+nskip
   mw=mw+1

endwhile

device,/close ; Z-buffer

if (plot_choice eq 'x') then begin
  set_plot,dname
  tv,img
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  filepng='pldt.png'
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.charsize=pcharsize

end
