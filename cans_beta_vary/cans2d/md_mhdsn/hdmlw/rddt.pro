narg=-1

dacgetparam,'params.txt','gm',gm
dacgetparam,'params.txt','ix',ix
dacgetparam,'params.txt','jx',jx
dacget0s,'t.dac',t,narg=narg
nx=n_elements(t)
dacget1d,'x.dac',x
dacget1d,'z.dac',z
dacget2s,'ro.dac',ro,narg=narg
dacget2s,'pr.dac',pr,narg=narg
dacget2s,'vx.dac',vx,narg=narg
dacget2s,'vz.dac',vz,narg=narg
dacget2s,'bx.dac',bx,narg=narg
dacget2s,'bz.dac',bz,narg=narg
dacget2s,'ay.dac',ay,narg=narg

ayr=ay
for i=0,ix-1 do begin
  ayr[i,*,*]=sqrt(abs(ay[i,*,*])*x[i])
endfor

te=pr/ro*gm

end
