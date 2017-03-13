narg=-1

dacgetparam,'params.txt','gm',gm
dacgetparam,'params.txt','ix',ix
dacgetparam,'params.txt','jx',jx
dacgetparam,'params.txt','kx',kx
dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,'x.dac',x
dacget1d,'y.dac',y
dacget1d,'z.dac',z
dacget3s,'ro.dac',ro,narg=narg
dacget3s,'pr.dac',pr,narg=narg
dacget3s,'vx.dac',vx,narg=narg
dacget3s,'vy.dac',vy,narg=narg
dacget3s,'vz.dac',vz,narg=narg
dacget3s,'bx.dac',bx,narg=narg
dacget3s,'by.dac',by,narg=narg
dacget3s,'bz.dac',bz,narg=narg

te=pr/ro*gm

end
