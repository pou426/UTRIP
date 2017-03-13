narg=-1

dacgetparam,'params.txt','gm',gm
dacgetparam,'params.txt','ix',ix
dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,'x.dac',x
dacget1d,'ro.dac',ro
dacget1d,'sc.dac',sc
dacget1s,'pr.dac',pr,narg=narg

te=pr
for n=0,nx-1 do te(*,n)=pr(*,n)/ro*gm

end
