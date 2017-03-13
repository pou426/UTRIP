G= 6.67d-8 ; cm^3 s^-2 g^-1
Msun=2.0d33   ; g
m_0=1.673d-24 ; g
n0_0=1.d8 ; cm^(-3)

rr0= 1.d11  ; cm
vv0= sqrt(G*Msun/rr0)  ; cm/sec
tt0= rr0/vv0  ; sec
ro0= m_0*n0_0 ; g/cm^3
pr0= ro0*vv0^2 ; erg/cm^3
bb0= sqrt(pr0) ; Gauss

zetac_0=1.104d11 ; cm
nu_0=2.42d20 ;cm^3 s^-2 g^(-1/3)
eta_0=5.24d-8 ; s^(-2)
a0_0=1.5d21 ; G cm^2
lambda_0=(5.54d11)^(-1) ; cm^(-1)

zetac = zetac_0/rr0
nu = nu_0/(vv0^2/ro0^(1/3.))
eta = eta_0/(1.d0/tt0^2)
a0 = a0_0/(sqrt(pr0)*rr0^2)
lambda = lambda_0/(1.d0/rr0)

print,'zetac=',zetac
print,'nu=',nu
print,'eta=',eta
print,'A0=',a0
print,'lambda=',lambda

end
