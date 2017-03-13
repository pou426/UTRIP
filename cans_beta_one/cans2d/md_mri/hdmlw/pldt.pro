
if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
iskip=2 & jskip=4 & scale=100 & limit=1.d-4
xrange=[min(x),max(x)] & yrange=[min(z),max(z)]
nskip=2

nlevels=10
dmin=0.5d-4 & dmax=2.d-4 & dlevel=(dmax-dmin)/nlevels
levels2=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
iwx=2 & jwx=3
if (plot_interactive eq 1) then $
read,' Plot columns & rows ? : ',iwx,jwx

ans='by'
;if (plot_interactive eq 1) then $
;read,' Variable for color-maps ? (by) : ',ans

nstart=10
if (plot_interactive eq 1) then $
read,'  start step ? : ',nstart

nlevels=32
color_index=(256/nlevels)*indgen(nlevels)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
case ans of
  'by' : begin
     data1=by
     dmin=-1.d-3 & dmax=1.d-3 & dlevel=(dmax-dmin)/nlevels
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
     ,levels1=levels1,clr_index=color_index $
     ,ay(*,*,n) $
     ,levels2=levels2 $
     ,vx(*,*,n),vz(*,*,n) $
     ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,sample=0 $
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
