narg=-1

dacgetparam,'params.txt','cs2',cs2
dacgetparam,'params.txt','csca2',csca2
dacgetparam,'params.txt','ix',ix
dacgetparam,'params.txt','jx',jx
dacgetparam,'params.txt','thini',thini
dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,'x.dac',x
dacget1d,'y.dac',y
dacget2s,'ro.dac',ro,narg=narg
dacget2s,'vx.dac',vx,narg=narg
dacget2s,'vy.dac',vy,narg=narg
dacget2s,'vz.dac',vz,narg=narg
dacget2s,'bx.dac',bx,narg=narg
dacget2s,'by.dac',by,narg=narg
dacget2s,'bz.dac',bz,narg=narg
dacget2s,'az.dac',az,narg=narg

end
