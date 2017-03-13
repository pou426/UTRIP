narg=-1

dacgetparam,'params.txt','gm',gm
dacgetparam,'params.txt','thini',thini
dacgetparam,'params.txt','phini',phini
dacgetparam,'params.txt','ix',ix
dacgetparam,'params.txt','jx',jx
dacgetparam,'params.txt','kx',kx
dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,'x.dac',x
dacget1d,'y.dac',y
dacget1d,'z.dac',z
dacget3d,'ro.dac',ro
dacget3s,'pr.dac',pr,narg=narg

te=pr
for n=0,nx-1 do te(*,*,*,n)=pr(*,*,*,n)/ro*gm

end
