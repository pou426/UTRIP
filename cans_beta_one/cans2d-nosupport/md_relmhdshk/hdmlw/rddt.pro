narg=-1

dacgetparam,'params.txt','mtest',mtest
dacgetparam,'params.txt','mdirection',mdirection
dacgetparam,'params.txt','gm',gm
dacgetparam,'params.txt','ix',ix
dacgetparam,'params.txt','jx',jx
dacgetparam,'params.txt','thini',thini
dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,'x.dac',x
dacget1d,'y.dac',y
dacget2s,'ro.dac',ro,narg=narg
dacget2s,'pr.dac',pr,narg=narg
dacget2s,'vx.dac',vx,narg=narg
dacget2s,'vy.dac',vy,narg=narg
dacget2s,'bx.dac',bx,narg=narg
dacget2s,'by.dac',by,narg=narg
dacget2s,'az.dac',az,narg=narg

gl=1./sqrt(1.-vx^2-vy^2)
te=pr/ro*gm

end
