
if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
iskip=10 & jskip=5 & scale=0.2 & limit=0.1
xrange=[min(x),max(x)] & yrange=[min(y),max(y)]
nskip=2

levels2=-60+4.*findgen(20)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
iwx=7 & jwx=1
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
     data1=ro
     levels1=0.2*findgen(16)
  end
  'pr' : begin
     data1=pr
     levels1=0.5*findgen(16)
  end
  'te' : begin
     data1=te
     dmin=0.d0 & dmax=4.d0 & dlevel=(dmax-dmin)/nlevels
     levels1=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)
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
     ,levels1=levels1,clr_index=clr_index $
     ,az(*,*,n) $
     ,levels2=levels2 $
     ,vx(*,*,n),vy(*,*,n) $
     ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,sample=0 $
     ,xcord=x,ycord=y
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
