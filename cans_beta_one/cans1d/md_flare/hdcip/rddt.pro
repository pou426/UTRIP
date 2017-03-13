narg=-1

dacgetparam,'params.txt','gm',gm
dacgetparam,'params.txt','rlnml',rlnml
dacgetparam,'params.txt','ronml',ronml
dacgetparam,'params.txt','tenml',tenml
dacgetparam,'params.txt','ix',ix
dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,'x.dac',x
dacget1d,'sc.dac',sc
dacget1s,'ro.dac',ro,narg=narg
dacget1s,'pr.dac',pr,narg=narg
dacget1s,'vx.dac',vx,narg=narg

te=pr/ro*gm

end
