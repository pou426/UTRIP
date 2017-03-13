!p.multi=[0,2,2]
!x.range=[1,4]

nt=nx-1
;nt=30
;nt=20

plot,x,ro(*,0),yr=[0.99,1.03],title='Density'
for n=0,nx-1 do oplot,x,ro(*,n)
;oplot,x,ro(*,nt)

plot,x,-vx(*,nt),title='Vx',/yl,/xl
for n=1,nx-1 do oplot,x,-vx(*,n)
;oplot,x,vx(*,nt)

plot,x,pr(*,0),yr=[0,0.002],title='Pressure'
for n=0,nx-1 do oplot,x,pr(*,n)
;oplot,x,pr(*,nt)

plot,x,gl(*,0),yr=[0.9,1.4],title='Lorentz Factor'
for n=0,nx-1 do oplot,x,gl(*,n)
;oplot,x,gl(*,nt)

!p.multi=0
end
