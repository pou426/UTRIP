narg=-1

dacgetparam,'params.txt','cs2',cs2
dacgetparam,'params.txt','ix',ix
dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,'x.dac',x
dacget1d,'bx.dac',bx
dacget1s,'ro.dac',ro,narg=narg
dacget1s,'vx.dac',vx,narg=narg
dacget1s,'vy.dac',vy,narg=narg
dacget1s,'by.dac',by,narg=narg


end
