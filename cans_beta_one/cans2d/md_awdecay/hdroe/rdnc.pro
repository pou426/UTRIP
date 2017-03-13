file='out.nc'

nnarray=-1

ncgetos,file,'cs2',cs2
ncgetos,file,'csca2',csca2
ncgetos,file,'thini',thini

ncgeto1,file,'x',x         & ix=n_elements(x)
ncgeto1,file,'y',y         & jx=n_elements(y)
ncgetss,file,nnarray,'t',t & nx=n_elements(t)

ncgets2,file,nnarray,'ro',ro
ncgets2,file,nnarray,'vx',vx
ncgets2,file,nnarray,'vy',vy
ncgets2,file,nnarray,'bx',bx
ncgets2,file,nnarray,'by',by
ncgets2,file,nnarray,'vz',vz
ncgets2,file,nnarray,'bz',bz
ncgets2,file,nnarray,'az',az


end
