iskip=5 & jskip=5 & scale=1 & limit=0.2
xrange=[0,3] & yrange=[0,3]
nskip=1

levels2=0.01*findgen(30)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
iwx=0 & jwx=0
read,' Plot columns & rows ? : ',iwx,jwx

ans=' '
read,' Variable for color-maps ? (ro,pr,te) : ',ans

nstart=0
read,'  start step ? : ',nstart

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
case ans of
  'ro' : begin
     data1=alog10(ro)
     levels1=-3+0.25*findgen(16)
  end
  'pr' : begin
     data1=pr
     levels1=0.003*findgen(16)
  end
  'te' : begin
     data1=te
     levels1=0.1*findgen(15)
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
     ,levels1=levels1,clr_index=16*bindgen(16) $
     ,ay(*,*,n) $
     ,levels2=levels2 $
     ,vx(*,*,n),vz(*,*,n) $
     ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,sample=0 $
     ,xcord=x,ycord=z
   put_time,t(n)

   n=n+nskip
   mw=mw+1

endwhile


!p.multi=0

end
