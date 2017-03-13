outdir='./moviedt'
if (not file_test(outdir,/directory)) then file_mkdir,outdir

iwx=3 & jwx=2

set_plot,'z'
xsize=iwx*240 & ysize=jwx*240
device,set_resolution=[xsize,ysize]
device,set_character_size=[6,9]

!p.charsize=2.0
!x.margin=[5,2]

nnarray=lindgen(nx)
xrange=[min(x),max(x)]

;;;

mx=5
namear=strarr(mx)
dminar=fltarr(mx)
dmaxar=fltarr(mx)
dlogar=fltarr(mx)

m=0 & namear[m]='vx' & dminar[m]=-6.00d0 & dmaxar[m]=6.00d0
m=1 & namear[m]='vy' & dminar[m]=-4.00d0 & dmaxar[m]=11.00d0
m=2 & namear[m]='ro' & dminar[m]=20.^(-9.0d0) & dmaxar[m]=1.3 & dlogar[m]=1
m=3 & namear[m]='pr' & dminar[m]=10.^(-7.0d0) & dmaxar[m]=10.^(-4.0d0) & dlogar[m]=1
m=4 & namear[m]='by' & dminar[m]=-5.00d-3 & dmaxar[m]=5.0d-3

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

    if (m eq 0) then put_time,t[n],charsize=1.5,xp=0.01,yp=0.97

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
