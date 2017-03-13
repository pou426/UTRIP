narg=-1
;narg=indgen(31)*10

dacgetparam,'params.txt','ca0',ca0
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
dacget3s,'ro.dac',ro,narg=narg
dacget3s,'vx.dac',vx,narg=narg
dacget3s,'vy.dac',vy,narg=narg
dacget3s,'vz.dac',vz,narg=narg
dacget3s,'bx.dac',bx,narg=narg
dacget3s,'by.dac',by,narg=narg
dacget3s,'bz.dac',bz,narg=narg

end
