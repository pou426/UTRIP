narg=-1
;narg=[15]
;narg=indgen(7)*10
dir='./'

mpest='.'+string(0,form='(i4.4)')
dacgetparam,dir+'params.txt'+mpest   ,'mpex',mpex
dacgetparam,dir+'params.txt'+mpest   ,'margin',margin
dacgetparam,dir+'params.txt'+mpest   ,'gm',gm
dacgetparam,dir+'params.txt'+mpest   ,'thini',thini
dacgetparam,dir+'params.txt'+mpest   ,'phini',phini
dacgetparam,dir+'params.txt'+mpest   ,'ebx',ebx
dacgetparam,dir+'params.txt'+mpest   ,'eby',eby
dacgetparam,dir+'params.txt'+mpest   ,'ebz',ebz
dacget0s,dir+'t.dac'+mpest,t,narg=narg
nx=n_elements(t)

mpestar='.'+string(indgen(mpex),form='(i4.4)')

ipear=intarr(mpex)
jpear=intarr(mpex)
kpear=intarr(mpex)
for mpe=0,mpex-1 do begin
  mpest=mpestar[mpe]
  dacgetparam,dir+'params.txt'+mpest   ,'ipe',ipe
  dacgetparam,dir+'params.txt'+mpest   ,'jpe',jpe
  dacgetparam,dir+'params.txt'+mpest   ,'kpe',kpe
  ipear[mpe]=ipe
  jpear[mpe]=jpe
  kpear[mpe]=kpe
endfor

files=dir+'x.dac'+mpestar[uniq(ipear,sort(ipear))]
dacget1d,files,x,margin=margin
files=dir+'y.dac'+mpestar[uniq(jpear,sort(jpear))]
dacget1d,files,y,margin=margin
files=dir+'z.dac'+mpestar[uniq(kpear,sort(kpear))]
dacget1d,files,z,margin=margin
ix=n_elements(x)
jx=n_elements(y)
kx=n_elements(z)

if 1 then begin
files=dir+'ro.dac'+mpestar
dacget4s,files,ro,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
files=dir+'pr.dac'+mpestar
dacget4s,files,pr,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
files=dir+'vx.dac'+mpestar
dacget4s,files,vx,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
files=dir+'vy.dac'+mpestar
dacget4s,files,vy,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
files=dir+'vz.dac'+mpestar
dacget4s,files,vz,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
files=dir+'bx.dac'+mpestar
dacget4s,files,bx,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
files=dir+'by.dac'+mpestar
dacget4s,files,by,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
files=dir+'bz.dac'+mpestar
dacget4s,files,bz,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
endif

te=pr/ro*gm
pm=(bx^2+by^2+bz^2)/!pi/8.e0

end
