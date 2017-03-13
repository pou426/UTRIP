narg=-1

dacgetparam,'params.txt','gm',gm
dacgetparam,'params.txt','ix',ix
dacgetparam,'params.txt','jx',jx
dacgetparam,'params.txt','kx',kx
dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,'r.dac',r
dacget1d,'ph.dac',ph
dacget1d,'z.dac',z
dacget3s,'ro.dac',ro,narg=narg
dacget3s,'pr.dac',pr,narg=narg
dacget3s,'vr.dac',vr,narg=narg
dacget3s,'vph.dac',vph,narg=narg
dacget3s,'vz.dac',vz,narg=narg
dacget3s,'br.dac',br,narg=narg
dacget3s,'bph.dac',bph,narg=narg
dacget3s,'bz.dac',bz,narg=narg

te=pr/ro*gm

end
