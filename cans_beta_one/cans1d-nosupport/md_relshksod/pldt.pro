!p.multi=[0,2,2]
!x.style=1
nt0=n_elements(t)-1

plot,x,de(*,0),title='De',yrange=[0,12]*1.d4,linest=1
oplot,x,de(*,nt0)


plot,x,vx(*,0),yr=[-0.5,1.5]*2.d-3,title='vx',linest=1
oplot,x,vx(*,nt0)

plot,x,eis(*,0),title='Specific Internal Energy',yrange=[0,2.2]*2.d-5,linest=1
oplot,x,eis(*,nt0)

!p.multi=0
!x.style=0

end
