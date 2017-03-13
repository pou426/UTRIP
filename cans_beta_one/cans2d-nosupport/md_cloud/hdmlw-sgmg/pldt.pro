if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

iskip=5 & jskip=5 & scale=0.05 & limit=1.d-4
xrange=[min(x),max(x)] & yrange=[min(y),max(y)]
nskip=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
iwx=3 & jwx=2
if (plot_interactive eq 1) then $
read,' Plot columns & rows ? : ',iwx,jwx

;ans=' '
;read,' Variable for color-maps ? (ro,pr,te) : ',ans
ans='ro'

nstart=0
if (plot_interactive eq 1) then $
read,'  start step ? : ',nstart

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
case ans of
  'ro' : begin
     data1=alog10(ro)
     levels1=-1+0.2*findgen(16)
  end
endcase

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

mwx=iwx*jwx

!p.multi=[0,iwx,jwx]

n=nstart
mw=0

while((mw lt mwx) and (n lt nx)) do begin

    fccnve $
     ,novector=0,nocontour=1 $
     ,xrange=xrange,yrange=yrange $
     ,data1(*,*,n) $
     ,levels1=levels1,clr_index=16*bindgen(16) $
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
