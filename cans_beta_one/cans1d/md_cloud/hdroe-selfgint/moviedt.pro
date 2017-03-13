outdir='./moviedt'
if (not file_test(outdir,/directory)) then file_mkdir,outdir

iwx=3 & jwx=1

set_plot,'z'
xsize=iwx*240 & ysize=jwx*240
device,set_resolution=[xsize,ysize]
device,set_character_size=[6,9]

!p.charsize=1.2
!x.margin=[5,2]

nnarray=lindgen(nx)
xrange=[min(x),max(x)]

;;;

mx=3
namear=strarr(mx)
dminar=fltarr(mx)
dmaxar=fltarr(mx)
dlogar=fltarr(mx)

m=0 & namear[m]='vx' & dminar[m]=-1.30d1 & dmaxar[m]=1.30d1
m=1 & namear[m]='ro' & dminar[m]=0.10d0 & dmaxar[m]=1.00d2 & dlogar[m]=1
m=2 & namear[m]='gx' & dminar[m]=-1.30d2 & dmaxar[m]=1.30d2

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

for m=0,mx-1 do begin

    name=namear[m]
    void=execute('data='+name+'[*,'+strtrim(n,2)+']')
    dmin=dminar[m] & dmax=dmaxar[m]

    plot,x,data ,/xst,xrange=xrange,yrange=[dmin,dmax],title=name,ylog=dlogar[m]
    oplot,!x.crange,[0,0],linestyle=1

    if (m eq 0) then put_time,t[n],charsize=1.5,xp=0.01,yp=0.95

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
