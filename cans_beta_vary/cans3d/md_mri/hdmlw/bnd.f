c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,vz,bx,by,bz
     &   ,dvy,dt,yrgn,y,ix,jx,kx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension y(jx),ysa(jx),ysb(jx)
      dimension js0(jx)
c----------------------------------------------------------------------|      
      call bdpery(margin,margin,pr,ix,jx,kx)
      call bdpery(margin,margin,ro,ix,jx,kx)
      call bdpery(margin,margin,vx,ix,jx,kx)
      call bdpery(margin,margin,vy,ix,jx,kx)
      call bdpery(margin,margin,vz,ix,jx,kx)
      call bdpery(margin,margin,bx,ix,jx,kx)
      call bdpery(margin,margin,by,ix,jx,kx)
      call bdpery(margin,margin,bz,ix,jx,kx)

      call bdperx(margin,margin,pr,ix,jx,kx)
      call bdperx(margin,margin,ro,ix,jx,kx)
      call bdperx(margin,margin,vx,ix,jx,kx)
      call bdperx(margin,margin,vy,ix,jx,kx)
      call bdperx(margin,margin,vz,ix,jx,kx)
      call bdperx(margin,margin,bx,ix,jx,kx)
      call bdperx(margin,margin,by,ix,jx,kx)
      call bdperx(margin,margin,bz,ix,jx,kx)

      if (.false.) then 
      do k=1,kx
      do j=1,jx
      do i=1,margin
        vy(i,j,k) =vy(i,j,k)-dvy
      enddo
      enddo
      enddo
      do k=1,kx
      do j=1,jx
      do i=1,margin
        vy(ix-margin+i,j,k) =vy(ix-margin+i,j,k)+dvy
      enddo
      enddo
      enddo

      else
c----------------------------------------------------------------------|      
c 'sliding block' boundary 
c     m=0 : i=1  side
c     m=1 : i=ix side

c     yrgn=y(jx-margin)-y(margin)
      ibnd0=1+margin
      ibnd1=ix-margin

      do m=0,1

      if (m.eq.0) then
        msign=-1
        ibnd00=ibnd0
        ibnd11=ibnd1+1
      else
        msign=+1
        ibnd00=ibnd1
        ibnd11=ibnd0-1
      endif


      dly=msign*dvy*dt

      if (dly.ge.0.d0) then
        dly0=mod(dly,yrgn)
      else
        dly0=mod(dly,yrgn)+yrgn
      endif

      do j=1,jx
        ysa(j)=y(j)+dly0
        ysb(j)=y(j)+dly0-yrgn
      enddo

      jab=1
      do j=1,jx
        if (y(j).gt.ysa(1)) then
          jab=j
          goto 100
        endif
      enddo
100   continue

      do j=1,jab-1
      js0(j)=jx
      do js=1,jx
        if (ysb(js).ge.y(j)) then
          js0(j)=js-1
          goto 210
        endif
      enddo
210   continue
      enddo
      do j=jab,jx
      js0(j)=jx
      do js=1,jx
        if (ysa(js).ge.y(j)) then
          js0(j)=js-1
          goto 220
        endif
      enddo
220   continue
      enddo

      do j=1,jx
        if (j.lt.jab) then
          dys=ysb(js0(j)+1)-ysb(js0(j))
          dy1=ysb(js0(j)+1)-y(j)
          dy0=y(j)-ysb(js0(j))
        else
          dys=ysa(js0(j)+1)-ysa(js0(j))
          dy1=ysa(js0(j)+1)-y(j)
          dy0=y(j)-ysa(js0(j))
        endif

        do i=1,margin
        do k=1,kx
          ro(ibnd00+msign*i,j,k)=(ro(ibnd11+msign*i,js0(j)+1,k)*dy0
     &                           +ro(ibnd11+msign*i,js0(j)  ,k)*dy1)/dys
          pr(ibnd00+msign*i,j,k)=(pr(ibnd11+msign*i,js0(j)+1,k)*dy0
     &                           +pr(ibnd11+msign*i,js0(j)  ,k)*dy1)/dys
          vx(ibnd00+msign*i,j,k)=(vx(ibnd11+msign*i,js0(j)+1,k)*dy0
     &                           +vx(ibnd11+msign*i,js0(j)  ,k)*dy1)/dys
          vy(ibnd00+msign*i,j,k)=(vy(ibnd11+msign*i,js0(j)+1,k)*dy0
     &                           +vy(ibnd11+msign*i,js0(j)  ,k)*dy1)/dys
     &                           +msign*dvy
          vz(ibnd00+msign*i,j,k)=(vz(ibnd11+msign*i,js0(j)+1,k)*dy0
     &                           +vz(ibnd11+msign*i,js0(j)  ,k)*dy1)/dys
          bx(ibnd00+msign*i,j,k)=(bx(ibnd11+msign*i,js0(j)+1,k)*dy0
     &                           +bx(ibnd11+msign*i,js0(j)  ,k)*dy1)/dys
          by(ibnd00+msign*i,j,k)=(by(ibnd11+msign*i,js0(j)+1,k)*dy0
     &                           +by(ibnd11+msign*i,js0(j)  ,k)*dy1)/dys
          bz(ibnd00+msign*i,j,k)=(bz(ibnd11+msign*i,js0(j)+1,k)*dy0
     &                           +bz(ibnd11+msign*i,js0(j)  ,k)*dy1)/dys
        enddo
        enddo

      enddo

      enddo
      endif


      call bdperz(margin,margin,pr,ix,jx,kx)
      call bdperz(margin,margin,ro,ix,jx,kx)
      call bdperz(margin,margin,vx,ix,jx,kx)
      call bdperz(margin,margin,vy,ix,jx,kx)
      call bdperz(margin,margin,vz,ix,jx,kx)
      call bdperz(margin,margin,bx,ix,jx,kx)
      call bdperz(margin,margin,by,ix,jx,kx)
      call bdperz(margin,margin,bz,ix,jx,kx)

      return
      end
