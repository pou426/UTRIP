Ca=1.
Cs=0.8
k = 1.

for i=1,360 do begin 
theta=i*!pi/180.+1.e-3
kz = k*cos(theta)
ky=  k*sin(theta)

;slow mode

omega=(0.5)^(0.5)*k*(Cs^2+Ca^2-((Cs^2+Ca^2)^2-4*Cs^2*Ca^2*(kz/k)^2)^0.5)^0.5
A= (2*omega^2-(Cs^2+Ca^2)*k^2)
vy=(ky/omega)*((Cs^2+Ca^2)*omega^2 - Ca^2*Cs^2*kz^2)/A
vz=(kz/omega)*((Cs^2+Ca^2)*omega^2 - Ca^2*Cs^2*(kz^2+k^2))/A

xyouts,vz*100+250,vy*100+260,'o',/dev

;fast mode

omega=(0.5)^(0.5)*k*(Cs^2+Ca^2+((Cs^2+Ca^2)^2-4*Cs^2*Ca^2*(kz/k)^2)^0.5)^0.5
A= (2*omega^2-(Cs^2+Ca^2)*k^2)
vy=(ky/omega)*((Cs^2+Ca^2)*omega^2 - Ca^2*Cs^2*kz^2)/A
vz=(kz/omega)*((Cs^2+Ca^2)*omega^2 - Ca^2*Cs^2*(kz^2+k^2))/A

xyouts,vz*100+250,vy*100+260,'x',/dev

endfor

xyouts,180,100, 'Ca=1, Cs=0.8',size=2.,/dev

end
