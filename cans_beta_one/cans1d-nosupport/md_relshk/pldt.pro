!p.multi=[0,2,2]
!x.style=1
nt0=n_elements(t)-1

plot,x,vx(*,0),yr=[0.0,0.8],title='vx',linest=1
oplot,!x.crange,[1,1]*0.69,linest=2
oplot,x,vx(*,nt0)

plot,x,eis(*,0),title='Specific Internal Energy',yrange=[0,2.2],linest=1
oplot,x,eis(*,nt0)

plot,x,de(*,0),title='De',yrange=[0,12],linest=1
oplot,!x.crange,[1,1]*3.55,linest=2
oplot,!x.crange,[1,1]*6.85,linest=2
oplot,x,de(*,nt0)

plot,x,pr(*,0),title='Pressure',yrange=[0,14],linest=1
oplot,!x.crange,[1,1]*1.384,linest=2
oplot,x,pr(*,nt0)

!p.multi=0
!x.style=0

end
