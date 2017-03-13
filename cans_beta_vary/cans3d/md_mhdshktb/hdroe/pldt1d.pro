
if (n_elements(plot_choice) eq 0) then plot_choice='x'

if (n_elements(plot_interactive) eq 0) then plot_interactive=1

if (plot_choice eq 'png') then begin
  set_plot,'z'
  xsize=640 & ysize=480
  device,set_resolution=[xsize,ysize]
  device,set_character_size=[6,9]
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
x0=-0.5*sin(thini)*cos(phini) & x1= 0.5*sin(thini)*cos(phini)
y0=-0.5*sin(thini)*sin(phini) & y1= 0.5*sin(thini)*sin(phini)
z0=-0.5*cos(thini)            & z1= 0.5*cos(thini)

xi=x0+(x1-x0)/22.*findgen(23)
yi=y0+(y1-y0)/22.*findgen(23)
zi=z0+(z1-z0)/22.*findgen(23)

s=sqrt((xi-x0)^2+(yi-y0)^2+(zi-z0)^2)
s=s-max(s)/2

ires=interpol(findgen(ix),x,xi)
jres=interpol(findgen(jx),y,yi)
kres=interpol(findgen(kx),z,zi)

n=0
ro_s=interpolate(ro(*,*,*,n),ires,jres,kres)
pr_s=interpolate(pr(*,*,*,n),ires,jres,kres)
vx_s=interpolate(vx(*,*,*,n),ires,jres,kres)
vy_s=interpolate(vy(*,*,*,n),ires,jres,kres)
vz_s=interpolate(vz(*,*,*,n),ires,jres,kres)
bx_s=interpolate(bx(*,*,*,n),ires,jres,kres)
by_s=interpolate(by(*,*,*,n),ires,jres,kres)
bz_s=interpolate(bz(*,*,*,n),ires,jres,kres)
vs_s=vx_s*sin(thini)*cos(phini) $
    +vy_s*sin(thini)*sin(phini) $
    +vz_s*cos(thini)
bs_s=bx_s*sin(thini)*cos(phini) $
    +by_s*sin(thini)*sin(phini) $
    +bz_s*cos(thini)
vt_s=vx_s*ebx +vy_s*eby +vz_s*ebz
bt_s=bx_s*ebx +by_s*eby +bz_s*ebz
te_s=gm*pr_s/ro_s

ro_s0=ro_s
pr_s0=pr_s
te_s0=te_s
vs_s0=vs_s
vt_s0=vt_s
bs_s0=bs_s
bt_s0=bt_s

n=3
ro_s=interpolate(ro(*,*,*,n),ires,jres,kres)
pr_s=interpolate(pr(*,*,*,n),ires,jres,kres)
vx_s=interpolate(vx(*,*,*,n),ires,jres,kres)
vy_s=interpolate(vy(*,*,*,n),ires,jres,kres)
vz_s=interpolate(vz(*,*,*,n),ires,jres,kres)
bx_s=interpolate(bx(*,*,*,n),ires,jres,kres)
by_s=interpolate(by(*,*,*,n),ires,jres,kres)
bz_s=interpolate(bz(*,*,*,n),ires,jres,kres)
vs_s=vx_s*sin(thini)*cos(phini) $
    +vy_s*sin(thini)*sin(phini) $
    +vz_s*cos(thini)
bs_s=bx_s*sin(thini)*cos(phini) $
    +by_s*sin(thini)*sin(phini) $
    +bz_s*cos(thini)
vt_s=vx_s*ebx +vy_s*eby +vz_s*ebz
bt_s=bx_s*ebx +by_s*eby +bz_s*ebz
te_s=gm*pr_s/ro_s

!p.multi=[0,3,2]
!x.style=1

plot,s,pr_s0,title='Pr',yrange=[0,1.2],linest=1
oplot,s,pr_s

plot,s,ro_s0,title='De',yrange=[0,1.2],linest=1
oplot,s,ro_s

plot,s,te_s0,yrange=[0.,5],title='Te',linest=1
oplot,s,te_s

plot,s,vs_s0,yr=[-0.5,1],title='Vs',linest=1
oplot,s,vs_s

plot,s,vt_s0,yr=[-1.7,0.2],title='Vt',linest=1,/yst
oplot,s,vt_s

plot,s,bt_s0/sqrt(!pi*4),yr=[-1.2,1.2],title='Bt/sqrt(4pi)',linest=1,/yst
oplot,s,bt_s/sqrt(!pi*4)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (plot_choice eq 'png') then begin
  filepng='pldt1d.png'
  img=tvrd()
  tvlct,red,green,blue,/get
  write_png,filepng,img,red,green,blue
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
!x.style=0

end
