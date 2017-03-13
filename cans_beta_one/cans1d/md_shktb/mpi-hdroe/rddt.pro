narg=-1

dacgetparam,'params.txt.0000','gm',gm
dacgetparam,'params.txt.0000','ix',ix
dacgetparam,'params.txt.0000','mpex',mpex
dacgetparam,'params.txt.0000','ipex',ipex
dacgetparam,'params.txt.0000','margin',margin
dacgetparam,'params.txt.0000','ro1',ro1
dacgetparam,'params.txt.0000','pr1',pr1

dacget0s,'t.dac.0000',t,narg=narg
nx=n_elements(t)
files='x.dac.'+string(indgen(ipex),form='(i4.4)')
dacget1d,files,x,margin=margin
files='ro.dac.'+string(indgen(mpex),form='(i4.4)')
dacget1s,files,ro,narg=narg,margin=margin
files='pr.dac.'+string(indgen(mpex),form='(i4.4)')
dacget1s,files,pr,narg=narg,margin=margin
files='vx.dac.'+string(indgen(mpex),form='(i4.4)')
dacget1s,files,vx,narg=narg,margin=margin

te=pr/ro*gm

end
