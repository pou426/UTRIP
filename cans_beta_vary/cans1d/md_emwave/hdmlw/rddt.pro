narg=-1

dacgetparam,'params.txt','ix',ix
dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,'x.dac',x
dacget1s,'ey.dac',ey,narg=narg
dacget1s,'ez.dac',ez,narg=narg
dacget1s,'by.dac',by,narg=narg
dacget1s,'bz.dac',bz,narg=narg
end
