iwx=0 & jwy=0
read,' window size ? : ',iwx,jwy

ans=' '
read,'  color    ? (ro,pr,te,bz) : ',ans

;ans1=' '
;read,'  contour ? (az) : ',ans1
ans1='az'

mmstart=0
read,'  start step ? : ',mmstart

xrange=[0,5.0] & yrange=[-5,5]
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

;levels2=0.04+0.02*findgen(16)
levels2=0.02+0.072*findgen(16)
data2=az_an

xplot=0.5*findgen(33)
yplot=0.+0.5*findgen(33)
iskip=5 & jskip=5
scale=1.0 & limit=0.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

iw=iwx*jwy
dx0=1./iwx
dy0=1./jwy

!p.multi=[0,iwx,jwy]

mm=mmstart
i=0

costh=(fltarr(ix)+1.)#cos(y)
sinth=(fltarr(ix)+1.)#sin(y)

while((i lt iw) and (mm lt nx)) do begin

    ; Magnetic Field lines
    if (ans1 eq 'az') then begin
      data20=reform(data2(*,*,mm))*(x#sin(y))
    endif else begin
      data20=reform(data2(*,*,mm))
    endelse

    datavx=reform(bx_an(*,*,mm))*sinth+reform(by_an(*,*,mm))*costh
    datavy=reform(bx_an(*,*,mm))*costh-reform(by_an(*,*,mm))*sinth

    ; Polar Coordinate Mapping
    ;
    xmin=max(x)*sin(y(0))
    dxx=0.05
    dd1=transpose(polar_surface(data1(*,*,mm),x,y,/grid,spacing=[dxx,dxx] $
      ,missing=min(data1(*,*,mm))))
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

    printf,-1,mm,format='(i3,$)'
    fccnve $
     ,novector=0,nocontour=0 $
     ,xrange=xrange,yrange=yrange $
     ,dd1 $
     ,levels1=levels1,clr_index=clr $
     ,dd2 $
     ,levels2=levels2 $
     ,ddvx,ddvy $
;    ,xplot=xplot,yplot=yplot,scale=scale,limit=limit $
     ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,sample=0 $
     ,xcord=xx,ycord=yy

	x0=dx0*(i mod iwx)
	y0=dy0*((i-(i mod iwx))/iwx)

;  circle,0,0,1,thick=4
s0=-!pi/2+!pi/100.*findgen(101)
polyfill,cos(s0),sin(s0),color=0

   put_time,t(mm)
mm=mm+6
i=i+1
endwhile

print,' '

; back to default

!x.range=[0.,0.]
!y.range=[0.,0.]
!x.margin=[10.,3.]
!y.margin=[4.,2.]
!p.charsize=1.
!p.thick=1.
!p.multi=0


delvar,data1,data2
end
