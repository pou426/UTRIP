!p.charsize=2

x0=min(x)*cos(thini) & x1= min(x)*cos(thini)
y0=min(y)*sin(thini) & y1= max(y)*sin(thini)

slice_1dfrom2d,ro1d,s,ro,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,pr1d,s,pr,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,te1d,s,te,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,vx1d,s,vx,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,vy1d,s,vy,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,bx1d,s,bx,[x0,y0],[x1,y1],x=x,y=y
slice_1dfrom2d,by1d,s,by,[x0,y0],[x1,y1],x=x,y=y

vs1d=vx1d*cos(thini)+vy1d*sin(thini)
vn1d=-vx1d*sin(thini)+vy1d*cos(thini)
bs1d=bx1d*cos(thini)+by1d*sin(thini)
bn1d=-bx1d*sin(thini)+by1d*cos(thini)


!p.multi=[0,3,2]
!x.style=1

n=nx-1

plot,s,pr1d(*,0),title='Pr',yrange=[0,1.2]*1.d-4,linest=1
oplot,s,pr1d(*,n)

plot,s,ro1d(*,0),title='De',yrange=[0,1.2],linest=1
oplot,s,ro1d(*,n)

plot,s,te1d(*,0) ,yrange=[0.,5]*1.d-4,title='Te',linest=1
oplot,s,te1d(*,n)

plot,s,vs1d(*,0),yr=[-0.4,0.8]*0.02,title='Vs',linest=1
oplot,s,vs1d(*,n)

plot,s,vn1d(*,0),yr=[-1.7,0.2]*0.01,title='Vn',linest=1,/yst
oplot,s,vn1d(*,n)

plot,s,bn1d(*,0)/sqrt(!pi*4),yr=[-1.2,1.2]*0.01,title='Bn/sqrt(4pi)',linest=1,/yst
oplot,s,bn1d(*,n)/sqrt(!pi*4)


!p.charsize=1
!p.multi=0
!x.style=0

end
