outdir='./moviedt'
if (not file_test(outdir,/directory)) then file_mkdir,outdir

iwx=4 & jwx=2

set_plot,'z'
xsize=iwx*240 & ysize=jwx*240
device,set_resolution=[xsize,ysize]
device,set_character_size=[6,9]

!p.charsize=1.2
!x.margin=[5,2]

iskip=10 & jskip=5 & scale=0.1 & limit=1.d-4
xrange=[min(x),max(x)] & yrange=[min(y),max(y)]
nlevels=32
color_index=(256/nlevels)*indgen(nlevels)


nnarray=lindgen(nx)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


mx=7
namear=strarr(mx)
dminar=fltarr(mx)
dmaxar=fltarr(mx)

m=0 & namear[m]='ro' & dminar[m]=0.0d0 & dmaxar[m]=8.0d0
m=1 & namear[m]='vx' & dminar[m]=-7.0d0 & dmaxar[m]=+7.0d0
m=2 & namear[m]='vy' & dminar[m]=-5.0d0 & dmaxar[m]=+5.0d0
m=3 & namear[m]='vz' & dminar[m]=-2.0d0 & dmaxar[m]=+2.0d0
m=4 & namear[m]='bx' & dminar[m]=-20.0d0 & dmaxar[m]=+20.0d0
m=5 & namear[m]='by' & dminar[m]=-25.0d0 & dmaxar[m]=+25.0d0
m=6 & namear[m]='bz' & dminar[m]=-8.0d0 & dmaxar[m]=+8.0d0

;for m=0,mx-1 do begin
;  name=namear[m]
;  void=execute('dminar['+strtrim(m,2)+']=min('+name+')')
;  void=execute('dmaxar['+strtrim(m,2)+']=max('+name+')')
;endfor

nmvx=n_elements(nnarray)

for nmv=0,nmvx-1 do begin

!p.multi=[0,iwx,jwx]

n=nnarray[nmv]

for m=0,mx-1 do begin
     name=namear[m]
     void=execute('data='+name+'[*,*,'+strtrim(n,2)+']')
     dmin=dminar(m) & dmax=dmaxar(m)
     dlevel=(dmax-dmin)/nlevels
     levels1=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)
    fccnve $
     ,novector=0,nocontour=1 $
     ,xrange=xrange,yrange=yrange $
     ,data $
     ,levels1=levels1,clr_index=color_index $
     ,vx[*,*,n],vy[*,*,n] $
     ,iskip=iskip,jskip=jskip,scale=scale,limit=limit $
     ,xcord=x,ycord=y,title=name
   put_time,t(nmv)
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
