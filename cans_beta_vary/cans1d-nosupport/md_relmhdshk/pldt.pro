bxar=bx#(dblarr(nx)+1)

ss=x
mx=ix

uss=gl*vx
prs=pr
ros=ro
vss=vx
vts=vy
bss=bxar
bts=by

theta=atan(vy,vx)
bparas=bxar*cos(theta)+by*sin(theta)
bperps=-bxar*sin(theta)+by*cos(theta)
bparacmv=bparas
bperpcmv=bperps/gl
pmcmv=(bperpcmv^2+bparacmv^2)/8./!pi

uss0=uss[0,0] & uss1=uss[mx-1,0]
prs0=prs[0,0] & prs1=prs[mx-1,0]
ros0=ros[0,0] & ros1=ros[mx-1,0]
bts0=bts[0,0] & bts1=bts[mx-1,0]
pmcmv0=pmcmv[0,0] & pmcmv1=pmcmv[mx-1,0]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

!x.style=1
!y.style=1

n=nx-1

case mtest of
1 : begin
  vw=0.5 & ussrange=[0.8,1.8] & prsrange=[0.,120.] & rosrange=[0.,5.]
end
2 : begin
  vw=0.0 & ussrange=[0.8,2.0] & prsrange=[0.,20.] & rosrange=[1.,3.]
end
3 : begin
  vw=0.2 & ussrange=[0.8,2.0] & prsrange=[0.,12.] & rosrange=[1.3,1.8]
end
4 : begin
  vw=0.0 & ussrange=[-2.5,0.] & prsrange=[0.,12.] & btsrange=[-1.,6.]
end
5 : begin
  vw=0.0 & ussrange=[-0.8,0.2] & prsrange=[0.,1.2] & btsrange=[-0.2,1.2]
end
6 : begin
  vssrange=[-0.004,0.008] & vtsrange=[-0.017,0.002]
  prsrange=[0.,1.2d-4] & rosrange=[0.,1.2]
  btsrange=[-0.012,0.012]
end
7 : begin
  vssrange=[-0.2,0.6] & vtsrange=[-0.6,0.2]
  prsrange=[0.,1.2] & rosrange=[0.,1.2]
  btsrange=[-1.2,1.2]
end
8 : begin
  ussrange=[-1.0,3.0] & prsrange=[0.3,2.d3] & rosrange=[0.,1.2]
end
9 : begin
  ussrange=[-0.5,2.0] & prsrange=[0.3,3.d2] & rosrange=[0.,1.2]
end
endcase

case 1 of
(mtest lt 4) : begin
xdc=vw*t[n]

!p.multi=[0,2,2]
plot,ss,uss[*,n],yr=ussrange,title='us',psym=4
oplot,[!x.crange[0],xdc],uss0*[1,1]
oplot,[xdc,!x.crange[1]],uss1*[1,1]
oplot,[xdc,xdc],[uss0,uss1]

plot,ss,prs[*,n],yr=prsrange,title='pr, pm_comv',psym=4
oplot,[!x.crange[0],xdc],prs0*[1,1]
oplot,[xdc,!x.crange[1]],prs1*[1,1]
oplot,[xdc,xdc],[prs0,prs1]

oplot,ss,pmcmv[*,n],psym=4
oplot,[!x.crange[0],xdc],pmcmv0*[1,1]
oplot,[xdc,!x.crange[1]],pmcmv1*[1,1]
oplot,[xdc,xdc],[pmcmv0,pmcmv1]

plot,ss,ros[*,n],yr=rosrange,title='ro',psym=4
oplot,[!x.crange[0],xdc],ros0*[1,1]
oplot,[xdc,!x.crange[1]],ros1*[1,1]
oplot,[xdc,xdc],[ros0,ros1]

end
(mtest ge 4 and mtest lt 6) : begin
xdc=vw*t[n]

!p.multi=[0,2,2]

plot,ss,uss[*,n],yr=ussrange,title='us',psym=4
oplot,[!x.crange[0],xdc],uss0*[1,1]
oplot,[xdc,!x.crange[1]],uss1*[1,1]
oplot,[xdc,xdc],[uss0,uss1]

plot,ss,prs[*,n],yr=prsrange,title='pr',psym=4
oplot,[!x.crange[0],xdc],prs0*[1,1]
oplot,[xdc,!x.crange[1]],prs1*[1,1]
oplot,[xdc,xdc],[prs0,prs1]

plot,ss,bts[*,n]/sqrt(!pi*4.),yr=btsrange,title='bt/sqrt(4*pi)',psym=4
oplot,[!x.crange[0],xdc],bts0/sqrt(!pi*4.)*[1,1]
oplot,[xdc,!x.crange[1]],bts1/sqrt(!pi*4.)*[1,1]
oplot,[xdc,xdc],[bts0,bts1]/sqrt(!pi*4.)

end
(mtest ge 6 and mtest lt 8) : begin
!p.multi=[0,2,3]

plot,ss,prs[*,n],yr=prsrange,title='pr',psym=4
plot,ss,ros[*,n],yr=rosrange,title='ro',psym=4
plot,ss,vss[*,n],yr=vssrange,title='vs',psym=4
plot,ss,vts[*,n],yr=vtsrange,title='vt',psym=4
plot,ss,bts[*,n]/sqrt(!pi*4.),yr=btsrange,title='bt/sqrt(4*pi)',psym=4

end
(mtest ge 8 and mtest lt 10) : begin
!p.multi=[0,2,2]

plot,ss,uss[*,n],yr=ussrange,title='ux',psym=4
plot,/ylog,ss,prs[*,n],yr=prsrange,title='pr',psym=4
oplot,ss,pmcmv[*,n]+prs[*,n],psym=4
plot,ss,ros[*,n],yr=rosrange,title='ro',psym=4

end
endcase


!p.multi=0
!x.style=0
!y.style=0

end
