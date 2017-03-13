narg=-1

dacgetparam,'params.txt','gm',gm
dacgetparam,'params.txt','zetac_0',zetac_0
dacgetparam,'params.txt','rnu_0',rnu_0
dacgetparam,'params.txt','eta_0',eta_0
dacgetparam,'params.txt','a0_0',a0_0
dacgetparam,'params.txt','rlambda_0',rlambda_0
dacgetparam,'params.txt','rMsun',rMsun
dacgetparam,'params.txt','rr0',rr0
dacgetparam,'params.txt','rn0_0',rn0_0

dacgetparam,'params.txt','ix',ix
dacgetparam,'params.txt','jx',jx

;dacget2d,'gx.dac',gx
;dacget2d,'gy.dac',gy

dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)

dacget1d,'x.dac',x
dacget1d,'y.dac',y
dacget2s,'ro.dac',ro,narg=narg
dacget2s,'pr.dac',pr,narg=narg
dacget2s,'vx.dac',vx,narg=narg
dacget2s,'vy.dac',vy,narg=narg
dacget2s,'vz.dac',vz,narg=narg
dacget2s,'bx.dac',bx,narg=narg
dacget2s,'by.dac',by,narg=narg
dacget2s,'bz.dac',bz,narg=narg
dacget2s,'az.dac',az,narg=narg

te=pr/ro*gm

end
