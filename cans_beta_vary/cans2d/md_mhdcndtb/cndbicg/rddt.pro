narg=-1

dacgetparam,'params.txt','gm',gm
dacgetparam,'params.txt','ix',ix
dacgetparam,'params.txt','jx',jx
dacgetparam,'params.txt','thini',thini
dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,'x.dac',x
dacget1d,'y.dac',y
dacget2d,'ro.dac',ro
dacget2d,'bx.dac',bx
dacget2d,'by.dac',by
dacget2d,'az.dac',az
dacget2s,'pr.dac',pr,narg=narg

te=pr
for n=0,nx-1 do te(*,*,n)=pr(*,*,n)/ro*gm

end
