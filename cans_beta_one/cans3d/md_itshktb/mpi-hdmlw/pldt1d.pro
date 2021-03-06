
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
vx_s=interpolate(vx(*,*,*,n),ires,jres,kres)
vy_s=interpolate(vy(*,*,*,n),ires,jres,kres)
vz_s=interpolate(vz(*,*,*,n),ires,jres,kres)
vs_s=vx_s*sin(thini)*cos(phini) $
    +vy_s*sin(thini)*sin(phini) $
    +vz_s*cos(thini)
ro_s0=ro_s & vs_s0=vs_s

n=3
ro_s=interpolate(ro(*,*,*,n),ires,jres,kres)
vx_s=interpolate(vx(*,*,*,n),ires,jres,kres)
vy_s=interpolate(vy(*,*,*,n),ires,jres,kres)
vz_s=interpolate(vz(*,*,*,n),ires,jres,kres)
vs_s=vx_s*sin(thini)*cos(phini) $
    +vy_s*sin(thini)*sin(phini) $
    +vz_s*cos(thini)


!p.multi=[0,2,1]

!x.style=1
nt0=n_elements(t)-1

plot,s,ro_s0,title='De',yrange=[0,1.2],linest=1
oplot,s,ro_s,psym=4

plot,s,vs_s0,yr=[-0.5,1.5],title='vx',linest=1
oplot,s,vs_s,psym=4

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
