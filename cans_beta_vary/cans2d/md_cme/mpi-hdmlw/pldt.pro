if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

iskip=5 & jskip=5 & scale=1.0 & limit=0.
xrange=[0,5.0] & yrange=[-5,5]
nskip=1

levels2=0.02+0.072*findgen(16)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

iwx=3 & jwx=2
if (plot_interactive eq 1) then $
read,' Plot columns & rows ? : ',iwx,jwx

ans='ro'
if (plot_interactive eq 1) then $
read,' Variable for color-maps ? (ro,pr,te,bz) : ',ans

nstart=0
if (plot_interactive eq 1) then $
read,'  start step ? : ',nstart

nlevels=32
color_index=(256/nlevels)*indgen(nlevels)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
case ans of
'ro' : begin
   data1=alog10(ro)
     dmin=-3.d0 & dmax=1.0d0 & dlevel=(dmax-dmin)/nlevels
     levels1=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)
end
'pr' : begin
   data1=pr
   levels1=0.0004*findgen(16)
end
'te' : begin
   data1=te
   levels1=0.21*findgen(16)
end
'bz' : begin
   data1=bz
   levels1=0.05*findgen(16)
end
endcase


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

mwx=iwx*jwx

!p.multi=[0,iwx,jwx]

n=nstart
mw=0

costh=(fltarr(ix)+1.)#cos(y)
sinth=(fltarr(ix)+1.)#sin(y)

while((mw lt mwx) and (n lt nx)) do begin

    ; Magnetic Field lines
      data20=reform(az(*,*,n))*(x#sin(y))

    datavx=reform(bx(*,*,n))*sinth+reform(by(*,*,n))*costh
    datavy=reform(bx(*,*,n))*costh-reform(by(*,*,n))*sinth

    ; Polar Coordinate Mapping
    ;
    xmin=max(x)*sin(y(0))
    dxx=0.05
    dd1=transpose(polar_surface(data1(*,*,n),x,y,/grid,spacing=[dxx,dxx] $
      ,missing=min(data1(*,*,n))))
    dd2=transpose(polar_surface(data20,x,y,/grid,spacing=[dxx,dxx] $
      ,missing=min(data20)))
    ddvx=transpose(polar_surface(datavx,x,y,/grid,spacing=[dxx,dxx] $
      ,missing=0.))
    ddvy=transpose(polar_surface(datavy,x,y,/grid,spacing=[dxx,dxx] $
      ,missing=0.))
    sz=size(dd1)
    xx=xmin+dxx*(findgen(sz(1)))
    yy=-max(x)+dxx*(findgen(sz(2)))

    yyar=(fltarr(sz(1))+1.)#yy
    xxar=xx#(fltarr(sz(2))+1.)
    rrar=sqrt(xxar^2+yyar^2)
    whr=where(rrar lt 1.)
    ddvx(whr)=0.
    ddvy(whr)=0.

    fccnve $
     ,novector=0,nocontour=0 $
     ,xrange=xrange,yrange=yrange $
     ,dd1 $
     ,levels1=levels1,clr_index=color_index $
     ,dd2 $
     ,levels2=levels2 $
     ,ddvx,ddvy $
     ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,sample=0 $
     ,xcord=xx,ycord=yy

   ;  circle,0,0,1,thick=4
   s0=-!pi/2+!pi/100.*findgen(101)
   polyfill,cos(s0),sin(s0),color=0

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
; back to default

!p.multi=0


end
