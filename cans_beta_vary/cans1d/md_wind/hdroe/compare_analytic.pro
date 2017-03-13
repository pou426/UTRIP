function wind_analytic_func,x,y

common wind_analytic_param,gm

gm=1.05

dydx = 4.*(5.-3*gm-(gm-1)*exp(y)+2*(2*gm-3.)/exp(x)) $
         /(-(5-3*gm)+(gm+1)*exp(y)-4*(gm-1)/exp(x))

return,dydx
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wind_analytic,ro,pr,vx,x,gm,rstar,mdot,kk,ee

common wind_analytic_param,gmc

gmc=gm

ix=301

xl=-2.+0.01*findgen(ix)
yl=fltarr(ix)

whr=where(xl ge 0.)
icrit=min(whr)

yl(icrit)=0.
dydx=(-4*(gm-1)+2*sqrt(2*(5-3*gm)))/(gm+1)

for i=icrit+1,ix-1 do begin

  h=xl(i)-xl(i-1)

  y0=rk4(yl(i-1),dydx,xl(i-1),h,'wind_analytic_func')
  yl(i)=y0

  dydx=wind_analytic_func(xl(i),yl(i))

endfor

dydx=(-4*(gm-1)+2*sqrt(2*(5-3*gm)))/(gm+1)
for i=icrit-1,0,-1 do begin

  h=-(xl(i+1)-xl(i))

  y0=rk4(yl(i+1),dydx,xl(i+1),h,'wind_analytic_func')
  yl(i)=y0

  dydx=wind_analytic_func(xl(i),yl(i))

endfor

x0=exp(xl)
vx0=sqrt(exp(yl))

rcrit=(5-3*gm)/4/(gm-1)/gm*rstar^2/ee
cscrit=sqrt(rstar^2/gm/2/rcrit)

x=rcrit*x0
vx=cscrit*vx0

ro=mdot/(4*!pi*vx*x^2)
pr=kk*ro^gm

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

i=100
n=nx-1
kk=pr[i,n]/ro[i,n]^gm
mdot=ro[i,n]*(4*!pi*vx[i,n]*x[i]^2)
ee=-1.d0/gm*rstar^2/x[i]+pr[i,n]/ro[i,n]*gm/(gm-1.d0)+0.5d0*vx[i,n]^2

wind_analytic,roan,pran,vxan,xan,gm,rstar,mdot,kk,ee

csan=sqrt(gm*pran/roan)

ro0=ro(*,n)
pr0=pr(*,n)
vx0=vx(*,n)
te0=te(*,n)

!p.multi=[0,2,2]

plot,xan,vxan,linest=1,title='Velocity & SoundSpeed'
oplot,xan,csan,linest=1
oplot,x,vx0,psym=4
oplot,x,sqrt(te0)

plot,xan,roan,linest=1,/yl,title='Density'
oplot,x,ro0,psym=4

plot,xan,pran,linest=1,/yl,title='Pressure'
oplot,x,pr0,psym=4

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

end
