      implicit real (a-h,o-z)
      parameter(ix=208)
      double precision x(ix),t,gm
      double precision ro(ix),pr(ix),vx(ix),te(ix)
      real xr(ix),data(ix)
c---------------------------------------------------------------------

      nd=4

c---------------------------------------------------------------------

      mstatus=nf_open('out.nc',0,idf0)
      call ncgeto1(idf0,'x',x,ix)
      call ncgetos(idf0,'gm',gm)
      mstatus=nf_inq_varid(idf0,'ro',idro0)
      call ncgets1(idf0,idro0,nd,ro,ix)
      mstatus=nf_inq_varid(idf0,'pr',idpr0)
      call ncgets1(idf0,idpr0,nd,pr,ix)
      mstatus=nf_inq_varid(idf0,'vx',idvx0)
      call ncgets1(idf0,idvx0,nd,vx,ix)

c---------------------------------------------------------------------

c     gm=1.4d0
c     mf_x=11
c     call dacopnr1d(mf_x,'x.dac',mtype,ix0)
c     read(mf_x) x
c     mf_t=10
c     call dacopnr0s(mf_t,'t.dac',mtype,nx0)
c     mf_ro=12
c     call dacopnr1s(mf_ro,'ro.dac',mtype,ix0,nx0)
c     mf_pr=13
c     call dacopnr1s(mf_pr,'pr.dac',mtype,ix0,nx0)
c     mf_vx=14
c     call dacopnr1s(mf_vx,'vx.dac',mtype,ix0,nx0)
c     do n=1,nd
c       read(mf_t) t
c       read(mf_ro) ro
c       read(mf_pr) pr
c       read(mf_vx) vx
c     enddo

c---------------------------------------------------------------------

      do i=1,ix
        te(i)=pr(i)/ro(i)*gm
      enddo

c---------------------------------------------------------------------

      mstatus=pgopen('/xwin')

      call pgsubp(2,2)
      call pgsch(1.8)

      do i=1,ix
        xr(i)=x(i)
        data(i)=pr(i)
      enddo
      CALL PGENV(-0.5, 0.5, 0., 1.2,  0,  0)
      CALL PGLINE(ix, xr, data)
      CALL pglab('x','Pr','')

      do i=1,ix
        xr(i)=x(i)
        data(i)=ro(i)
      enddo
      CALL PGENV(-0.5, 0.5, 0., 1.2,  0,  0)
      CALL PGLINE(ix, xr, data)
      CALL pglab('x','Ro','')

      do i=1,ix
        xr(i)=x(i)
        data(i)=te(i)
      enddo
      CALL PGENV(-0.5, 0.5, 0., 2.0,  0,  0)
      CALL PGLINE(ix, xr, data)
      CALL pglab('x','Te','')

      do i=1,ix
        xr(i)=x(i)
        data(i)=vx(i)
      enddo
      CALL PGENV(-0.5, 0.5, -0.5, 1.5,  0,  0)
      CALL PGLINE(ix, xr, data)
      CALL pglab('x','Vx','')

      CALL PGCLOS

      END
