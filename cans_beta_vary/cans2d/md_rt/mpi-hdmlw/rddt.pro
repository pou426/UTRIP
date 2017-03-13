narg=-1
;narg=[15]
;narg=indgen(7)*10
dir='./'

mpest='.'+string(0,form='(i4.4)')
dacgetparam,dir+'params.txt'+mpest   ,'mpex',mpex
dacgetparam,dir+'params.txt'+mpest   ,'margin',margin
dacgetparam,dir+'params.txt'+mpest   ,'gm',gm
dacgetparam,dir+'params.txt'+mpest   ,'thini',thini
dacget0s,dir+'t.dac'+mpest,t,narg=narg
nx=n_elements(t)

mpestar='.'+string(indgen(mpex),form='(i4.4)')

ipear=intarr(mpex)
jpear=intarr(mpex)
for mpe=0,mpex-1 do begin
  mpest=mpestar[mpe]
  dacgetparam,dir+'params.txt'+mpest   ,'ipe',ipe
  dacgetparam,dir+'params.txt'+mpest   ,'jpe',jpe
  ipear[mpe]=ipe
  jpear[mpe]=jpe
endfor

files=dir+'x.dac'+mpestar[uniq(ipear,sort(ipear))]
dacget1d,files,x,margin=margin
files=dir+'y.dac'+mpestar[uniq(jpear,sort(jpear))]
dacget1d,files,y,margin=margin
ix=n_elements(x)
jx=n_elements(y)

if 1 then begin
files=dir+'ro.dac'+mpestar
dacget4s,files,ro,narg=narg,margin=margin,ipear=ipear,jpear=jpear
files=dir+'pr.dac'+mpestar
dacget4s,files,pr,narg=narg,margin=margin,ipear=ipear,jpear=jpear
files=dir+'vx.dac'+mpestar
dacget4s,files,vx,narg=narg,margin=margin,ipear=ipear,jpear=jpear
files=dir+'vy.dac'+mpestar
dacget4s,files,vy,narg=narg,margin=margin,ipear=ipear,jpear=jpear
endif

ro=reform(ro)
pr=reform(pr)
vx=reform(vx)
vy=reform(vy)
te=pr/ro*gm

end
