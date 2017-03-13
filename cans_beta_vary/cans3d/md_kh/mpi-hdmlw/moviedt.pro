outdir='./moviedt'
if (not file_test(outdir,/directory)) then file_mkdir,outdir

iwx=3 & jwx=2

dname=!d.name

xdw=240 & ydw=240
xsize=iwx*xdw & ysize=jwx*ydw

img=bytarr(xsize,ysize)

pcharsize=!p.charsize & !p.charsize=2.0
xmargin=!x.margin & !x.margin=[5,2]

nnarray=lindgen(nx)
xrange=[min(x),max(x)] & yrange=[min(y),max(y)] & zrange=[min(z),max(z)]

ax=30 & az=30
position=[0.01,0.01,0.99,0.99,0.26,0.74]
position=[0.01,0.01,0.99,0.99,0.01,0.99]

;;;
nlevels=32
color_index=(256/nlevels)*indgen(nlevels)
  
  xslice=0
  yslice=0.
  zslice=0

mx=5
namear=strarr(mx)
dminar=fltarr(mx)
dmaxar=fltarr(mx)
dlogar=fltarr(mx)

m=0 & namear[m]='vx' & dminar[m]=-0.50d0 & dmaxar[m]=0.50d0
m=1 & namear[m]='vy' & dminar[m]=-1.00d0 & dmaxar[m]=1.00d0
m=2 & namear[m]='vz' & dminar[m]=-0.50d0 & dmaxar[m]=0.50d0
m=3 & namear[m]='ro' & dminar[m]=0.76d0 & dmaxar[m]=1.04d0
m=4 & namear[m]='pr' & dminar[m]=0.54d0 & dmaxar[m]=0.64d0

;for m=0,mx-1 do begin
;  name=namear[m]
;  void=execute('dminar['+strtrim(m,2)+']=min('+name+')')
;  void=execute('dmaxar['+strtrim(m,2)+']=max('+name+')')
;endfor
;;;

;;;

iskip=4 & jskip=8 & scale=0.3 & limit=1.d-1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

nmvx=n_elements(nnarray)

for nmv=0,nmvx-1 do begin
!p.multi=[0,iwx,jwx]

n=nnarray[nmv]
mw=0

for m=0,mx-1 do begin

set_plot,'z'
device,set_resolution=[xdw,ydw]
device,set_character_size=[6,9]

    name=namear[m]
    if (dlogar[m] eq 1) then begin
      void=execute('data=alog10('+name+'[*,*,*,'+strtrim(n,2)+'])')
    endif else begin
      void=execute('data='+name+'[*,*,*,'+strtrim(n,2)+']')
    endelse
    dmin=dminar(m) & dmax=dmaxar(m)
    dlevel=(dmax-dmin)/nlevels
    levels1=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; plot

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

  if 0 then begin
  ; B-lines
  mfx=6 & radius=3
  s0=2.*!pi/mx*findgen(mx)
  y0=0. &  z0=-10.
  zplot=z0+radius*cos(s0)
  yplot=radius*sin(s0)

  xplot=[fltarr(mx)-30.,fltarr(mx)+30.]
  yplot=[yplot,yplot]
  zplot=[zplot,zplot]

  intline,bx[*,*,*,n],by[*,*,*,n],bz[*,*,*,n] $
    ,x,y,z,xplot,yplot,zplot $
    ,index,xline,yline,zline $
    ,n_points=3000
  mfx=n_elements(index)
  for mf=0,mfx-1 do begin
    m1=index(mf)
    if (mf eq 0) then m0=0 else m0=index(mf-1)+1
;   plots,/t3d,xline[m0:m1],yline[m0:m1],zline[m0:m1],/data
  if (m1 gt m0) then $
    plot_tube,/t3d,xline[m0:m1],yline[m0:m1],zline[m0:m1],radius=0.5
  endfor
  endif

  ; Y slice
  t3d,/yzexch
  zvalue=yslice*!y.s[1]+!y.s[0]
  zvalue=position[3]+0.01
  contour,/noerase $
    ,/t3d,position=position[[0,4,2,5]],zvalue=zvalue $
    ,reform(data[*,x2i(y,yslice),*]),x,z  $
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
    ,transpose(data[x2i(x,xslice),*,*]),z,y  $
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
    ,data[*,*,x2i(z,zslice)],x,y  $
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

  if (m eq 0) then put_time,t[n],charsize=1.5,xp=0.1,yp=0.9
  xyouts,/normal,0.1,0.1,name

  d=tvrd()
device,/close ; Z-buffer

  jw=mw/iwx
  iw=mw-jw*iwx
  img[iw*xdw:(iw+1)*xdw-1,(jwx-1-jw)*ydw:(jwx-jw)*ydw-1]=d

   mw=mw+1



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

endfor



filepng=outdir+'/'+string(nmv,format='(i3.3)')+'.png'
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue


endfor ;nmv

cd,current=current_dir
str=strmid(current_dir,strpos(current_dir,'cans',/reverse_search))
files=file_basename(file_search(outdir+'/*.png'))
movieindex,outdir+'/movie.html',files,title=str

!p.multi=0

end
!p.charsize=pcharsize
!x.margin=xmargin
set_plot,dname

end
