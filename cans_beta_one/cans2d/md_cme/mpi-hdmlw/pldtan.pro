iskip=5 & jskip=5 & scale=1.0 & limit=0.
xrange=[0,5.0] & yrange=[-5,5]
nskip=6

levels2=0.02+0.072*findgen(16)
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
   data1=alog10(ro_an)
   levels1=-3+0.25*findgen(16)
end
'pr' : begin
   data1=pr_an
   levels1=0.0004*findgen(16)
end
'te' : begin
   data1=te_an
   levels1=0.21*findgen(15)
end
'bz' : begin
   data1=bz_an
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
      data20=reform(az_an(*,*,n))*(x#sin(y))

    datavx=reform(bx_an(*,*,n))*sinth+reform(by_an(*,*,n))*costh
    datavy=reform(bx_an(*,*,n))*costh-reform(by_an(*,*,n))*sinth

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

    yyar=(fltarr(ix)+1.)#yy
    xxar=xx#(fltarr(jx)+1.)
    rrar=sqrt(xxar^2+yyar^2)
    whr=where(rrar lt 1.)
    ddvx(whr)=0.
    ddvy(whr)=0.

    fccnve $
     ,novector=0,nocontour=0 $
     ,xrange=xrange,yrange=yrange $
     ,dd1 $
     ,levels1=levels1,clr_index=clr $
     ,dd2 $
     ,levels2=levels2 $
     ,ddvx,ddvy $
     ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,sample=0 $
     ,xcord=xx,ycord=yy

   s0=-!pi/2+!pi/100.*findgen(101)
   polyfill,cos(s0),sin(s0),color=0

   put_time,t(n)

   n=n+nskip
   mw=mw+1
endwhile


; back to default

!p.multi=0


end
