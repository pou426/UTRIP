narg=-1

dacgetparam,'params.txt','mtest',mtest
dacgetparam,'params.txt','gm',gm
dacgetparam,'params.txt','ix',ix
dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,'x.dac',x
dacget1d,'bx.dac',bx
dacget1s,'ro.dac',ro,narg=narg
dacget1s,'pr.dac',pr,narg=narg
dacget1s,'vx.dac',vx,narg=narg
dacget1s,'vy.dac',vy,narg=narg
dacget1s,'by.dac',by,narg=narg

gl=1./sqrt(1.-vx^2-vy^2)
de=ro*gl
eis=pr/ro/(gm-1)
te=pr/ro*gm

end
