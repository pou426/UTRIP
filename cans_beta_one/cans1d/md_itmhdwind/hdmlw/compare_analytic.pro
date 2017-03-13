function newtfunc, in

common param1, gm,theta,omega

out=fltarr(6)

beta=in[0]
ee=in[1]

xx=in[2]  ; x-slow
yy=in[3]  ; y-slow

out[0]=beta/(2.d0*xx^4*yy^2)+theta/(gm-1.d0)*yy^(gm-1.d0)-1.d0/xx $
   +omega/2.d0*( (xx-1.d0/xx)^2/(yy-1.d0)^2 -xx^2)-ee
out[1]=-2.d0*beta/(xx^5*yy^2)+1.d0/xx^2-omega*xx $
     +omega*(xx-1.d0/xx)*(1.d0+1.d0/xx^2)/(yy-1.d0)^2
out[2]=-beta/(xx^4*yy^3)+theta*yy^(gm-2.d0) $
     -omega*(xx-1.d0/xx)^2/(yy-1.d0)^3

xx=in[4]  ; x-fast
yy=in[5]  ; y-fast

out[3]=beta/(2.d0*xx^4*yy^2)+theta/(gm-1.d0)*yy^(gm-1.d0)-1.d0/xx $
   +omega/2.d0*( (xx-1.d0/xx)^2/(yy-1.d0)^2 -xx^2)-ee
out[4]=-2.d0*beta/(xx^5*yy^2)+1.d0/xx^2-omega*xx $
     +omega*(xx-1.d0/xx)*(1.d0+1.d0/xx^2)/(yy-1.d0)^2
out[5]=-beta/(xx^4*yy^3)+theta*yy^(gm-2.d0) $
     -omega*(xx-1.d0/xx)^2/(yy-1.d0)^3

return,out
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function newtfun1, in

common param1, gm,theta,omega
common param2, beta,ee,xx

out=fltarr(1)

yy=in[0]

out[0]=beta/(2.d0*xx^4*yy^2)+theta/(gm-1.d0)*yy^(gm-1.d0)-1.d0/xx $
   +omega/2.d0*( (xx-1.d0/xx)^2/(yy-1.d0)^2 -xx^2)-ee

return,out
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; main routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

common param1, gm,theta,omega
common param2, beta,ee,xx

if (n_elements(plot_choice) eq 0) then plot_choice='x'
if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  filepng='compare_analytic.png'
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
gm=1.2d0

omegastart=0.01d0
thetastart=0.21d0
mox=241
mtx=201

domega=(0.25d0-omegastart)/(mox-1.d0)
dtheta=(0.50d0-thetastart)/(mtx-1.d0)

thetaar=thetastart+dtheta*dindgen(mtx)
omegaar=omegastart+domega*dindgen(mox)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

theta=thetaar[0]
omega=omegaar[0]

betag=2.d0^(2.d0*(3.d0-2*gm)/(gm-1.d0))*(2.d0*(gm-1.d0))^(-3.d0)* $
    ( (theta/(gm-1.d0)-1.d0) / (5.d0-3.d0*gm) )^( (5.d0-3.d0*gm)/(gm-1.d0) )
eeg=theta/(gm-1.d0)-1.d0
xsg=1.d0-4.d0*omega*(gm-1.d0)^2.d0/(1.d0-omega)^3.d0
ysg=1.d0+4.d0*omega*(gm-1.d0)/(1.d0-omega)^2.d0
xfg=(2.d0^(gm+1.d0)*theta^2.d0*betag^(gm-1.d0))^(-1.d0/(5.d0-3.d0*gm))
yfg=(16.d0*theta^3.d0*betag)^(1.d0/(5.d0-3.d0*gm))
in=[betag,eeg,xsg,ysg,xfg,yfg]

result=newton(in,'newtfunc',double=1,check=check)
in=result
;out=newtfunc(in)
;print,in[2],in[4]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for mo=1,mox-1 do begin
  omega=omegaar[mo]
  result=newton(in,'newtfunc',double=1,check=check)
  in=result
; y=newtfunc(in)
;  print,omega,in[2],in[4]
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for mt=1,mtx-1 do begin
  theta=thetaar[mt]
  result=newton(in,'newtfunc',double=1,check=check)
  in=result
  y=newtfunc(in)
;  print,theta,in[2],in[4]
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

beta=result[0]
ee=result[1]
xs=result[2]
ys=result[3]
xf=result[4]
yf=result[5]


xar0=(dindgen(301)+1.d0)*0.01d0
yar0=(dindgen(301)+1.1d0)*0.01d0
xar=xar0#(dblarr(301)+1.d0)
yar=(dblarr(301)+1.d0)#yar0

zz=beta/(2.*xar^4*yar^2)+theta/(gm-1.)*yar^(gm-1)-1./xar $
   +omega/2.*( (xar-1./xar)^2/(yar-1)^2 -xar^2)-ee

levels=[-0.1,0.,0.15,0.3,0.3848,0.5]
contour,zz,xar0,yar0,levels=levels,/xst,/yst
contour,zz,xar0,yar0,levels=[0],thick=3,/xst,/yst,/noerase

mx=231
x=[0.7+(dindgen(29)+1.0d0)*0.01d0,1.00d0+(dindgen(mx-29)+1.d0)*0.01d0]
y=dblarr(mx)

in=dblarr(1)
result=dblarr(1)

result[0]=3.0
for m=0,mx-1 do begin
  case 1 of
   (x[m] lt xs)                     : in[0]=result[0]
   ((x[m] ge xs) and (x[m] le 1.d0)): in[0]=1.01
   ((x[m] gt 1.d0) and (x[m] lt xf)): in[0]=0.99
   (x[m] ge xf)                     : in[0]=0.1
  endcase

  xx=x[m]
  result=newton(in,'newtfun1',double=1,check=check)
  y[m]=result[0]
  in[0]=result[0]
; print,x[m],y[m]
endfor

oplot,x,y,psym=4,symsize=3

;vy=(1./x-x*y)/(1.-y)
;plot,x,vy

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
