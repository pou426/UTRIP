c     ------------------------------------------------------------
      subroutine mglin(u,nx,ncycle,ng,memlen)
      integer nx,ncycle,npre,npost,ng,memlen
      double precision u(-1:nx,-1:nx)
c     this subroutine is given by `mutigrid methods for boundary value
c     problems. i. w. h.press & s.a.teukolsky, computer in physics, 5
c     514-519(1991)'
c     
c     variables:
c        ng = the number of grid's hierarchy.
c        npre = the number of pre-smoothing iteration.
c        npre = the number of post-smoothing iteration.
c     
c     uses addint,copy,fill0,interp,melloc,relax,resid,rstrct,slvsml
c     
      integer j,jcycle,jj,jpost,jpre,mem,nf,ngrid,nn,ires(ng),
     $     iro(ng),irhs(ng),iu(ng)
      double precision z(memlen)
      npre=2
      npost=2
      mem=0
      nn=nx/2
      ngrid=ng-1
      call melloc((nn+2)**2,iro,ng,ngrid,z,memlen,mem)
      call rstrct(z(iro(ngrid)),u,nn)
 1    if (nn.gt.2) then
         nn=nn/2
         ngrid=ngrid-1
         call melloc((nn+2)**2,iro,ng,ngrid,z,memlen,mem)
         call rstrct(z(iro(ngrid)),z(iro(ngrid+1)),nn)
         goto 1
      endif
      nn=2
      call melloc((nn+2)**2,iu,ng,1,z,memlen,mem)
      call melloc((nn+2)**2,irhs,ng,1,z,memlen,mem)
      call slvsml(z(iu(1)),z(iro(1)))
      ngrid=ng
      do j=2,ngrid
         nn=2*nn
         call melloc((nn+2)**2,iu,ng,j,z,memlen,mem)
         call melloc((nn+2)**2,irhs,ng,j,z,memlen,mem)
         call melloc((nn+2)**2,ires,ng,j,z,memlen,mem)
         call interp(z(iu(j)),z(iu(j-1)),nn)
         if (j.ne.ngrid) then
            call copy(z(irhs(j)),z(iro(j)),nn)
         else
            call copy(z(irhs(j)),u,nn)
         endif
         do jcycle=1,ncycle
            nf=nn
            do jj=j,2,-1
               do jpre=1,npre
                  call relax(z(iu(jj)),z(irhs(jj)),nf)
               enddo
               call resid(z(ires(jj)),z(iu(jj)),z(irhs(jj)),nf)
               nf=nf/2
               call rstrct(z(irhs(jj-1)),z(ires(jj)),nf)
               call fill0(z(iu(jj-1)),nf)
            enddo
            call slvsml(z(iu(1)),z(irhs(1)))
            nf=2
            do jj=2,j
               nf=2*nf
               call addint(z(iu(jj)),z(iu(jj-1)),z(ires(jj)),nf)
               do jpost=1,npost
                  call relax(z(iu(jj)),z(irhs(jj)),nf)
               enddo
            enddo
         enddo
      enddo
      call copy(u,z(iu(ngrid)),nx)
      return
      end
