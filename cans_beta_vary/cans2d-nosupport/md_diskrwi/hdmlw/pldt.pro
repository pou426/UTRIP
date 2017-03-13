iskip=5 & jskip=5 & scale=0.2 & limit=0.
xrange=[-2.0,2.0] & yrange=[-2,2]
nskip=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

iwx=0 & jwx=0
read,' Plot columns & rows ? : ',iwx,jwx

ans=' '
read,' Variable for color-maps ? (ro,pr,te,bz) : ',ans

nstart=0
read,'  start step ? : ',nstart

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
case ans of
'ro' : begin
   data1=ro
   levels1=1.+(0.01*findgen(16)-0.08)
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

    datavx=reform(vx(*,*,n))*sinth+reform(vy(*,*,n))*costh
    datavy=reform(vx(*,*,n))*costh-reform(vy(*,*,n))*sinth

    ; Polar Coordinate Mapping
    ;
    dxx=0.05
    dd1=transpose(polar_surface(data1(*,*,n),x,y,/grid,spacing=[dxx,dxx] $
      ,missing=min(data1(*,*,n))))
    ddvx=transpose(polar_surface(datavx,x,y,/grid,spacing=[dxx,dxx] $
      ,missing=0.))
    ddvy=transpose(polar_surface(datavy,x,y,/grid,spacing=[dxx,dxx] $
      ,missing=0.))
    sz=size(dd1)
    xx=-max(x)+dxx*(findgen(sz(1)))
    yy=-max(x)+dxx*(findgen(sz(2)))

    yyar=(fltarr(sz(1))+1.)#yy
    xxar=xx#(fltarr(sz(2))+1.)
    rrar=sqrt(xxar^2+yyar^2)
    whr=where(rrar lt 0.4)
    ddvx[whr]=0.
    ddvy[whr]=0.

    fccnve $
     ,novector=0,nocontour=1 $
     ,xrange=xrange,yrange=yrange $
     ,dd1 $
     ,levels1=levels1,clr_index=clr $
     ,ddvx,ddvy $
     ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,sample=0 $
     ,xcord=xx,ycord=yy

   ;  circle,0,0,1,thick=4
   s0=-!pi/2+!pi/100.*findgen(201)
   polyfill,0.4*cos(s0),0.4*sin(s0),color=0

   put_time,t(n)

   n=n+nskip
   mw=mw+1
endwhile


; back to default

!p.multi=0


end
