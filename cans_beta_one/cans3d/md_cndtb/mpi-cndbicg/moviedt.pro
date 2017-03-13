outdir='./moviedt'
if (not file_test(outdir,/directory)) then file_mkdir,outdir

iwx=1 & jwx=1

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

mx=1
namear=strarr(mx)
dminar=fltarr(mx)
dmaxar=fltarr(mx)

m=0 & namear[m]='te' & dminar[m]= 0.00d0 & dmaxar[m]=8.00d-1

;for m=0,mx-1 do begin
;  name=namear[m]
;  void=execute('dminar['+strtrim(m,2)+']=min('+name+')')
;  void=execute('dmaxar['+strtrim(m,2)+']=max('+name+')')
;endfor
;;;

;;;



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

    void=execute('data='+name+'[*,*,*,'+strtrim(n,2)+']')
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


  ; Y slice
  t3d,/yzexch
  zvalue=yslice*!y.s[1]+!y.s[0]
  zvalue=position[3]+0.01
  contour,/noerase $
    ,/t3d,position=position[[0,4,2,5]],zvalue=zvalue $
    ,reform(data[*,x2i(y,yslice),*]),x,z  $
    ,levels=levels1,xrange=xrange,yrange=zrange $
    ,xst=5,yst=5,/cell_fill,c_colors=color_index
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
  t3d,/xzexch

  ; Z slice
  zvalue=zslice*!z.s[1]+!z.s[0]
  zvalue=position[4]-0.01
  contour,/noerase $
    ,/t3d,position=position[[0,1,2,3]],zvalue=zvalue $
    ,data[*,*,x2i(z,zslice)],x,y  $
    ,levels=levels1,xrange=xrange,yrange=yrange $
    ,xst=5,yst=5,/cell_fill,c_colors=color_index

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
