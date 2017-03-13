narg=-1
dir='./'

dacgetparam,dir+'params.txt','gm',gm
dacgetparam,dir+'params.txt','ix',ix
dacget0s,dir+'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,dir+'x.dac',x
dacget1d,dir+'bx.dac',bx
dacget1s,dir+'ro.dac',ro,narg=narg
dacget1s,dir+'pr.dac',pr,narg=narg
dacget1s,dir+'vx.dac',vx,narg=narg
dacget1s,dir+'vy.dac',vy,narg=narg
dacget1s,dir+'by.dac',by,narg=narg

gl=1./sqrt(1.-vx^2-vy^2)
de=ro*gl
eis=pr/ro/(gm-1)
te=pr/ro*gm

end
