outdir='./moviedt'
if (not file_test(outdir,/directory)) then file_mkdir,outdir

iwx=4 & jwx=2

set_plot,'z'
xsize=iwx*150 & ysize=jwx*240
device,set_resolution=[xsize,ysize]
device,set_character_size=[6,9]

!p.charsize=2.0
!x.margin=[5,2]

nnarray=lindgen(nx)
;xrange=[min(x),max(x)] & yrange=[min(y),max(y)]
xrange=[0,5.0] & yrange=[-5,5]

;;;
nlevels=32
color_index=(256/nlevels)*indgen(nlevels)
  
mx=8
namear=strarr(mx)
dminar=fltarr(mx)
dmaxar=fltarr(mx)
dlogar=intarr(mx)

m=0 & namear[m]='ro' & dminar[m]=-3.00d0 & dmaxar[m]=0.00d0
m=1 & namear[m]='vx' & dminar[m]=-2.30d0 & dmaxar[m]=2.30d0
m=2 & namear[m]='vy' & dminar[m]=-0.10d0 & dmaxar[m]=0.10d0
m=3 & namear[m]='vz' & dminar[m]=-3.00d0 & dmaxar[m]=3.00d0
m=4 & namear[m]='pr' & dminar[m]=-4.00d0 & dmaxar[m]=-1.00d0
m=5 & namear[m]='bx' & dminar[m]=-0.30d0 & dmaxar[m]=0.30d0
m=6 & namear[m]='by' & dminar[m]=-0.20d0 & dmaxar[m]=0.20d0
m=7 & namear[m]='bz' & dminar[m]=-1.70d0 & dmaxar[m]=1.70d0
dlogar[0]=1
dlogar[4]=1

;for m=0,mx-1 do begin
;  name=namear[m]
;  void=execute('dminar['+strtrim(m,2)+']=min('+name+')')
;  void=execute('dmaxar['+strtrim(m,2)+']=max('+name+')')
;endfor
;;;

nlevels2=16
;dmax=max(az) & dmin=min(az)
;dlevel=(dmax-dmin)/nlevels2
;levels2=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels2-1)
levels2=0.02+0.072*findgen(nlevels2)
;;;

iskip=16 & jskip=16 & scale=1.0 & limit=1.d-2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

nmvx=n_elements(nnarray)

costh=(fltarr(ix)+1.)#cos(y)
sinth=(fltarr(ix)+1.)#sin(y)

for nmv=0,nmvx-1 do begin
!p.multi=[0,iwx,jwx]

n=nnarray[nmv]

for m=0,mx-1 do begin

    ; filled
    name=namear[m]
    if (dlogar[m] eq 1) then begin
      void=execute('data1=alog10('+name+'[*,*,'+strtrim(n,2)+'])')
    endif else begin
      void=execute('data1='+name+'[*,*,'+strtrim(n,2)+']')
    endelse
    dmin=dminar(m) & dmax=dmaxar(m)
    dlevel=(dmax-dmin)/nlevels
    levels1=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)

    ; Contour
    data2=reform(az[*,*,n])*(x#sin(y))

    ; Arrows
    datavx=reform(vx[*,*,n])*sinth+reform(vy[*,*,n])*costh
    datavy=reform(vx[*,*,n])*costh-reform(vy[*,*,n])*sinth

    ; Polar Coordinate Mapping
    ;
    xmin=max(x)*sin(y(0))
    dxx=0.05
    dd1=transpose(polar_surface(data1,x,y,/grid,spacing=[dxx,dxx] $
      ,missing=min(data1)))
    dd2=transpose(polar_surface(data2,x,y,/grid,spacing=[dxx,dxx] $
      ,missing=min(data2)))
    ddvx=transpose(polar_surface(datavx,x,y,/grid,spacing=[dxx,dxx] $
      ,missing=0.))
    ddvy=transpose(polar_surface(datavy,x,y,/grid,spacing=[dxx,dxx] $
      ,missing=0.))
    sz=size(dd1)
    xx=xmin+dxx*(findgen(sz[1]))
    yy=-max(x)+dxx*(findgen(sz[2]))

    yyar=(fltarr(sz[1])+1.)#yy
    xxar=xx#(fltarr(sz[2])+1.)
    rrar=sqrt(xxar^2+yyar^2)
    whr=where(rrar lt 1.)
    ddvx(whr)=0.
    ddvy(whr)=0.

    fccnve,xrange=xrange,yrange=yrange,xcord=xx,ycord=yy,title=name $
     ,nofill=0,dd1,levels1=levels1,clr_index=color_index $
     ,nocontour=0,dd2,levels2=levels2 $
     ,novector=0,ddvx,ddvy,iskip=iskip,jskip=jskip,scale=scale,limit=limit
    put_time,t(nmv)

   ;  circle,0,0,1,thick=4
   s0=-!pi/2+!pi/100.*findgen(101)
   polyfill,cos(s0),sin(s0),color=0

endfor


filepng=outdir+'/'+string(nmv,format='(i3.3)')+'.png'
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue


endfor ;nmv

cd,current=current_dir
str=strmid(current_dir,strpos(current_dir,'cans',/reverse_search))
files=file_basename(file_search(outdir+'/*.png'))
movieindex,outdir+'/movie.html',files,title=str

!p.multi=0

end
