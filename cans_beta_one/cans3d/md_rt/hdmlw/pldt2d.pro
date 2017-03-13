
if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xplot=0.5*findgen(33) & yplot=0.+0.5*findgen(33)
iskip=4 & jskip=10 & scale=0.1 & limit=1.d-1

mmskip=2
j=jx/2

xrange=[min(x),max(x)] & yrange=[min(z),max(z)]
nlevels=32
color_index=(256/nlevels)*indgen(nlevels)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
iwx=3 & jwy=2
if (plot_interactive eq 1) then $
read,' window size ? : ',iwx,jwy

ans='ro'
if (plot_interactive eq 1) then $
read,'  color    ? (ro,pr,te,vx,cy) : ',ans

;ans1=' '
;read,'  contour ? (az) : ',ans1
ans1='az'

mmstart=0
if (plot_interactive eq 1) then $
read,'  start step ? : ',mmstart
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
case ans of
'em' : begin
   data1=em
   levels1=0.02*findgen(16)
end
'ro' : begin
   data1=ro
   dmin=0.d0 & dmax=6.4d0 
   dlevel=(dmax-dmin)/nlevels
   levels1=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)
end
'pr' : begin
   data1=pr
   levels1=45*findgen(16)
end
'te' : begin
   data1=te
   levels1=6.*findgen(15)
end
've' : begin
   data1=ve
   levels1=0.1*findgen(15)
end
'vx' : begin
   data1=vx
   levels1=-0.69+0.1*findgen(15)
end
'vy' : begin
   data1=vy
   levels1=-0.69+0.1*findgen(15)
end
'cy' : begin
   data1=smooth(smooth(cy,2),2)
   levels1=-35.+5.*findgen(15)
   levels1=[-14+2*findgen(7),2+2*findgen(7)]
end
'et' : begin
   data1=et
   levels1=0.06*findgen(15)
end
'xr' : begin
;  data1=xr
;  levels1=0.25*findgen(15)
   data1=alog10(xr)
   levels1=0.2*findgen(15)
end
endcase

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

mwx=iwx*jwy

!p.multi=[0,iwx,jwy]

mm=mmstart
mw=0

while((mw lt mwx) and (mm lt nx)) do begin
    printf,-1,mm,format='(i3,$)'
    fccnve $
     ,novector=0,nocontour=1 $
     ,xrange=xrange,yrange=yrange $
     ,reform(data1(*,j,*,mm)) $
     ,levels1=levels1,clr_index=color_index $
;    ,data2(*,*,mm) $
;    ,levels2=levels2 $
     ,reform(vx(*,j,*,mm)),reform(vz(*,j,*,mm)) $
;    ,xplot=xplot,yplot=yplot,scale=scale,limit=limit $
     ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,sample=0 $
     ,xcord=x,ycord=z

   put_time,t(mm)
mm=mm+mmskip
mw=mw+1
endwhile

print,' '

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  filepng='pldt2d.png'
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; back to default

!x.range=[0.,0.]
!y.range=[0.,0.]
!x.margin=[10.,3.]
!y.margin=[4.,2.]
!p.charsize=1.
!p.thick=1.
!p.multi=0

end
