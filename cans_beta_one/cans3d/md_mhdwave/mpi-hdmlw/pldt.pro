if (n_elements(plot_choice) eq 0) then plot_choice='x'
if (n_elements(plot_interactive) eq 0) then plot_interactive=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
iskip=3 & jskip=3 & scale=100. & limit=1.d-4
xrange=[min(x),max(x)] & yrange=[min(y),max(y)] & zrange=[min(z),max(z)]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

iwx=3 & jwx=2
if (plot_interactive eq 1) then $
read,' Plot columns & rows ? : ',iwx,jwx

ans='ro'
if (plot_interactive eq 1) then $
read,' Variable for color-maps ? (ro,pr,te) : ',ans

nstart=0
if (plot_interactive eq 1) then $
read,'  start step ? : ',nstart

nskip=2
if (plot_interactive eq 1) then $
read,'  step interval ? : ',nskip

  xslice=0
  yslice=0.
  zslice=0

nlevels=32
color_index=(256/nlevels)*indgen(nlevels)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
case ans of
  'ro' : begin
     data1=ro
     dmin=0.9995d0 & dmax=1.0005d0
   end
endcase

dlevel=(dmax-dmin)/nlevels
levels1=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pcharsize=!p.charsize
!p.charsize=2.0

; plot

ax=30 & az=30
position=[0.01,0.01,0.99,0.99,0.26,0.74]
position=[0.01,0.01,0.99,0.99,0.01,0.99]


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


  ; Y slice
  t3d,/yzexch
  zvalue=yslice*!y.s[1]+!y.s[0]
  zvalue=position[3]+0.01
  contour,/noerase $
    ,/t3d,position=position[[0,4,2,5]],zvalue=zvalue $
    ,reform(data1[*,x2i(y,yslice),*,n]),x,z  $
    ,levels=levels1,xrange=xrange,yrange=zrange $
    ,xst=5,yst=5,/cell_fill,c_colors=color_index
  fccnve,/noerase,/noaxis,/nolabels $
    ,/t3d,position=position[[0,4,2,5]],zvalue=zvalue-0.01 $
    ,nofilled=1,novector=0,nocontour=1 $
    ,xcord=x,ycord=z $
    ,xrange=xrange,yrange=zrange $
    ,reform(vx[*,x2i(y,yslice),*,n]) $
    ,reform(vz[*,x2i(y,yslice),*,n]) $
    ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,v_thick=1.5,v_color=255
  t3d,/yzexch

  ; X slice
  t3d,/xzexch
  zvalue=xslice*!x.s[1]+!x.s[0]
  zvalue=position[2]+0.01
  contour,/noerase $
    ,/t3d,position=position[[4,1,5,3]],zvalue=zvalue $
    ,transpose(data1[x2i(x,xslice),*,*,n]),z,y  $
    ,levels=levels1,xrange=zrange,yrange=yrange $
    ,xst=5,yst=5,/cell_fill,c_colors=color_index
  fccnve,/noerase,/noaxis,/nolabels $
    ,/t3d,position=position[[4,1,5,3]],zvalue=zvalue-0.01 $
    ,nofilled=1,novector=0,nocontour=1 $
    ,xcord=z,ycord=y $
    ,xrange=zrange,yrange=yrange $
    ,transpose(reform(vz(x2i(x,xslice),*,*,n))) $
    ,transpose(reform(vy(x2i(x,xslice),*,*,n))) $
    ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,v_thick=1.5,v_color=255
  t3d,/xzexch

  ; Z slice
  zvalue=zslice*!z.s[1]+!z.s[0]
  zvalue=position[4]-0.01
  contour,/noerase $
    ,/t3d,position=position[[0,1,2,3]],zvalue=zvalue $
    ,data1[*,*,x2i(z,zslice),n],x,y  $
    ,levels=levels1,xrange=xrange,yrange=yrange $
    ,xst=5,yst=5,/cell_fill,c_colors=color_index
  fccnve,/noerase,/noaxis,/nolabels $
    ,/t3d,position=position[[0,1,2,3]],zvalue=zvalue+0.01 $
    ,nofilled=1,novector=0,nocontour=1 $
    ,xcord=x,ycord=y $
    ,xrange=xrange,yrange=yrange $
    ,vx[*,*,x2i(z,zslice),n] $
    ,vy[*,*,x2i(z,zslice),n] $
    ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,v_thick=1.5,v_color=255

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
