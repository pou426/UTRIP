outdir='./moviedt'
if (not file_test(outdir,/directory)) then file_mkdir,outdir

iwx=2 & jwx=2

set_plot,'z'
xsize=iwx*240 & ysize=jwx*240
device,set_resolution=[xsize,ysize]
device,set_character_size=[6,9]

!p.charsize=1.2
!x.margin=[5,2]

nnarray=lindgen(nx)
xrange=[min(x),max(x)]

;;;

mx=4
namear=strarr(mx)
dminar=fltarr(mx)
dmaxar=fltarr(mx)
dlogar=fltarr(mx)

m=0 & namear[m]='vx' & dminar[m]= 0.00d0 & dmaxar[m]=5.00d0
m=1 & namear[m]='vy' & dminar[m]= 0.00d0 & dmaxar[m]=5.00d-1
m=2 & namear[m]='ro' & dminar[m]=10.^(-7.0d0) & dmaxar[m]=1.00d0 & dlogar[m]=1
m=3 & namear[m]='by' & dminar[m]=-0.50d0 & dmaxar[m]=0.00d0

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

    plot,x,data ,/xst,xrange=xrange,yrange=[dmin,dmax],title=name,ylog=dlogar[m],xtitle='x'
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
