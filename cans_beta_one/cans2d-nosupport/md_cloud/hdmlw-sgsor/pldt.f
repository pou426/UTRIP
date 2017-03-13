      implicit real (a-h,o-z)
      parameter(ix=67,jx=67,ndx=4)
      double precision x(ix),y(jx),t
      double precision ro(ix,jx),vx(ix,jx),vy(ix,jx)
      real xr(ix),yr(jx),data(ix,jx),data1(ix,jx)
      parameter(ivx=20,jvx=20)
      real tr(6),trv(6)
      real datax(ivx,jvx),datay(ivx,jvx)
      integer nread(ndx)
      character tc*4

      nread(1)=1
      nread(2)=2
      nread(3)=4
      nread(4)=6

c---------------------------------------------------------------------
      mstatus=nf_open('out.nc',0,idf0)
      call ncgeto1(idf0,'x',x,ix)
      call ncgeto1(idf0,'y',y,jx)
      mstatus=nf_inq_varid(idf0,'ro',idro0)
      mstatus=nf_inq_varid(idf0,'vx',idvx0)
      mstatus=nf_inq_varid(idf0,'vy',idvy0)
      mstatus=nf_inq_varid(idf0,'t',idt0)
c---------------------------------------------------------------------

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
c     do i=1,ix
c       te(i)=pr(i)/ro(i)*gm
c     enddo

c---------------------------------------------------------------------

      tr(1)=x(1)
      tr(2)=x(2)-x(1)
      tr(3)=0.e0
      tr(4)=y(1)
      tr(5)=0.e0
      tr(6)=y(2)-y(1)

      trv(1)=x(1)
      trv(2)=(x(2)-x(1))*5
      trv(3)=0.e0
      trv(4)=y(1)
      trv(5)=0.e0
      trv(6)=(y(2)-y(1))*5

      mstatus=pgopen('/xwin')

      call pgsubp(2,2)

      do nd=1,ndx

      n=nread(nd)
      call ncgetss(idf0,idt0,n,t)
      call ncgets2(idf0,idro0,n,ro,ix,jx)
      call ncgets2(idf0,idvx0,n,vx,ix,jx)
      call ncgets2(idf0,idvy0,n,vy,ix,jx)

      do j=1,jx
      do i=1,ix
        data1(i,j)=log10(ro(i,j))
      enddo
      enddo

      vvlim=0.3

      do jv=1,jvx
      do iv=1,ivx
        vx0=vx(iv*5,jv*5)
        vy0=vy(iv*5,jv*5)
        if ((vx0**2+vy0**2).lt.vvlim) then
          datax(iv,jv)=0.e0
          datay(iv,jv)=0.e0
        else
          datax(iv,jv)=vx0
          datay(iv,jv)=vy0
        endif
      enddo
      enddo

      call pgsch(1.8)
      CALL PGENV( 0., 6.28, 0., 6.28,  0,  0)
      CALL PALETT(2, 1.0, 0.5)

      call pgimag(data1,ix,jx,1,ix,1,jx,-1.,2.,tr)
      call pgsch(0.5)
      call pgvect(datax,datay,ivx,jvx
     &         ,1,ivx,1,jvx,0.0,0,trv,0.e0)
      call pgsch(1.8)
      write(tc,'(f4.2)') t
      CALL pglab('x','y','t= '//tc)

      enddo

      CALL PGCLOS

      END

      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
C-----------------------------------------------------------------------
C Set a "palette" of colors in the range of color indices used by
C PGIMAG.
C-----------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
C
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
C
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
C
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
C
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C
      IF (TYPE.EQ.1) THEN
C        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
C        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
C        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
C        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
C        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
      END
