cd,'../cndsor'
.r rddt

ts=t
xs=x
ys=y
zs=z
;ros=ro
prs=pr
;vxs=vx
;vys=vy
;vzs=vz
;bxs=bx
;bys=by
;bzs=bz
;azs=az

cd,'../mpi-cndsor'
.r rddt

;dro=ro-ros & print,min(dro),max(dro)
dpr=pr-prs & print,min(dpr),max(dpr)
;dvx=vx-vxs & print,min(dvx),max(dvx)
;dvy=vy-vys & print,min(dvy),max(dvy)
;dvz=vz-vzs & print,min(dvz),max(dvz)
;dbx=bx-bxs & print,min(dbx),max(dbx)
;dby=by-bys & print,min(dby),max(dby)
;dbz=bz-bzs & print,min(dbz),max(dbz)
;daz=az-azs & print,min(daz),max(daz)

;end
