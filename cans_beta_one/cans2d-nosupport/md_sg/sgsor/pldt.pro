levels1=0.9+0.0125*findgen(17)
xrange=[0,2*!pi]
yrange=[0,2*!pi]
limit=0.
scale=3.0
iskip=5 & jskip=5


fccnve $
  ,xrange=xrange,yrange=yrange $
  ,ro,levels1=levels1,clr_index=clr $
  ,ro,levels2=levels1 $
  ,gx,gy $
  ,iskip=iskip,jskip=jskip,scale=scale,limit=limit,sample=0 $
  ,xcord=x,ycord=y

end
