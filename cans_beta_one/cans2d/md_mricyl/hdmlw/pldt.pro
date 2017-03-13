
if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
iskip=5 & jskip=5 & scale=0.1 & limit=0.2
;xplot=0.5*findgen(33) & yplot=0.+0.5*findgen(33)
xrange=[min(x),max(x)] & yrange=[min(z),max(z)]
nskip=4

nlevels=30
dmin=0. & dmax=0.3 & dlevel=(dmax-dmin)/nlevels
levels2=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
iwx=3 & jwy=2
if (plot_interactive eq 1) then $
read,' window size ? : ',iwx,jwy

ans='ro'
if (plot_interactive eq 1) then $
read,'  color    ? (ro,pr,te,vx,cy) : ',ans

nstart=0
if (plot_interactive eq 1) then $
read,'  start step ? : ',nstart

nlevels=32
color_index=(256/nlevels)*indgen(nlevels)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
case ans of
'em' : begin
   data1=em
   levels1=0.02*findgen(16)
end
'ro' : begin
   data1=alog10(ro)
     dmin=-3.d0 & dmax=1.d0 & dlevel=(dmax-dmin)/nlevels
     levels1=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)
end
'pr' : begin
   data1=pr
   levels1=0.0004*findgen(16)
end
'te' : begin
   data1=te
   levels1=0.21*findgen(15)
end
've' : begin
   data1=ve
   levels1=0.1*findgen(15)
end
'vx' : begin
   data1=vx
   levels1=-0.69+0.1*findgen(15)
end
'vz' : begin
   data1=vz
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
   data1=alog10(xr)
   levels1=0.2*findgen(15)
end
endcase

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

mwx=iwx*jwy

!p.multi=[0,iwx,jwy]

n=nstart
mw=0

while((mw lt mwx) and (n lt nx)) do begin

    fccnve $
     ,novector=0,nocontour=0 $
     ,xrange=xrange,yrange=yrange $
     ,data1(*,*,n) $
     ,levels1=levels1,clr_index=color_index $
     ,ay[*,*,n] $
     ,levels2=levels2 $
     ,vx(*,*,n),vz(*,*,n) $
;    ,xplot=xplot,yplot=yplot,scale=scale,limit=limit $
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
