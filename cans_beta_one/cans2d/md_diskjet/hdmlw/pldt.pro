if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
iskip=5 & jskip=5 & scale=0.2 & limit=0.2
xrange=[0,3] & yrange=[0,3]
xrange=[0,1.3] & yrange=[0,1.3]
nskip=2

levels2=0.04*findgen(30)
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

nlevels=32
color_index=(256/nlevels)*indgen(nlevels)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
case ans of
  'ro' : begin
     data1=alog10(ro)
     dmin=-3.d0 & dmax=1.d0 & dlevel=(dmax-dmin)/nlevels
     levels1=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)
  end
  'pr' : begin
     data1=pr
     levels1=0.003*findgen(16)
  end
  'te' : begin
     data1=alog10(te)
     levels1=-1.+0.1*findgen(16)
  end
endcase

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

mwx=iwx*jwx

!p.multi=[0,iwx,jwx]

n=nstart
mw=0

while((mw lt mwx) and (n lt nx)) do begin

    fccnve $
     ,novector=0,nocontour=0 $
     ,xrange=xrange,yrange=yrange $
     ,data1(*,*,n) $
;    ,levels1=levels1,clr_index=16*bindgen(16) $
     ,levels1=levels1,clr_index=color_index $
     ,ayr(*,*,n) $
     ,levels2=levels2 $
     ,vx(*,*,n),vz(*,*,n) $
     ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,sample=1,index_size=1 $
     ,xcord=x,ycord=z
   put_time,t(n)

   n=n+nskip
   mw=mw+1

endwhile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  filepng='pldt.png'
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

!p.multi=0

end
