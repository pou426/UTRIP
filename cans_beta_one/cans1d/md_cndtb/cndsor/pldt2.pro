n=30
dx=mkdelta(x)
q=total(te[*,n]*dx)

idxcnd=5/2.
;idxcnd=0.
ro0=1.
rkap0=0.1
kappa0=rkap0*(gm-1)/ro0

cndtb_analytic,f,xi,idxcnd

dxi=mkdelta(xi)
fttl=total(f*dxi)

xan=xi*(t[n]*kappa0*q^idxcnd)^(1./(idxcnd+2))
tean=(q^2/kappa0/t[n])^(1/(idxcnd+2.))*f/fttl
;dxan=mkdelta(xan)
;qan=total(tean*dxan)

!x.range=[0,1]

plot,x,te[*,n]
oplot,xan,tean,linest=2

!x.range=0

end
