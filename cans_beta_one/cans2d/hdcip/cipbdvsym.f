c======================================================================|
      subroutine cipbdvsym(vxm,vxdxm,vxdym
     &                    ,vym,vydxm,vydym
     &                    ,margin,ix,jx,mdir,mbnd,msign_n,msign_t)
c======================================================================|
c----------------------------------------------------------------------|
      implicit none
      integer margin,ix,jx
      integer mdir,mbnd,msign_n,msign_t
      integer marginn,margint
      double precision vxm,vxdxm,vxdym
      double precision vym,vydxm,vydym
      dimension vxm(ix,jx),vym(ix,jx)
      dimension vxdxm(ix,jx),vxdym(ix,jx)
      dimension vydxm(ix,jx),vydym(ix,jx)
c----------------------------------------------------------------------|

c     -- bottom
      if (mbnd.eq.0) then 
        marginn=margin-1
        margint=margin
c     -- top
      else
        marginn=margin
        margint=margin
      endif

c     -- x-direction
      if (mdir.eq.0) then 

c       -- Vn(bnd_out) = +Vn(bnd_in)
        if (msign_n.eq.+1) then 
          call bdsmpx(mbnd,marginn,vxm,ix,jx)
          call bdsmnx(mbnd,marginn,vxdxm,ix,jx)
          call bdsmpx(mbnd,marginn,vxdym,ix,jx)
c       -- Vn(bnd_out) = -Vn(bnd_in)
        else
          call bdsmnx(mbnd,marginn,vxm,ix,jx)
          call bdsmpx(mbnd,marginn,vxdxm,ix,jx)
          call bdsmnx(mbnd,marginn,vxdym,ix,jx)
        endif

c       -- Vt(bnd_out) = +Vt(bnd_in)
        if (msign_t.eq.+1) then 
          call bdsppx(mbnd,margint,vym,ix,jx)
          call bdspnx(mbnd,margint,vydxm,ix,jx)
          call bdsppx(mbnd,margint,vydym,ix,jx)
c       -- Vt(bnd_out) = -Vt(bnd_in)
        else
          call bdspnx(mbnd,margint,vym,ix,jx)
          call bdsppx(mbnd,margint,vydxm,ix,jx)
          call bdspnx(mbnd,margint,vydym,ix,jx)
        endif

      endif
c     -- x-direction

c     -- y-direction
      if (mdir.eq.1) then 

c       -- Vn(bnd_out) = +Vn(bnd_in)
        if (msign_n.eq.+1) then 
          call bdsmpy(mbnd,marginn,vym,ix,jx)
          call bdsmny(mbnd,marginn,vydym,ix,jx)
          call bdsmpy(mbnd,marginn,vydxm,ix,jx)
c       -- Vn(bnd_out) = -Vn(bnd_in)
        else
          call bdsmny(mbnd,marginn,vym,ix,jx)
          call bdsmpy(mbnd,marginn,vydym,ix,jx)
          call bdsmny(mbnd,marginn,vydxm,ix,jx)
        endif

c       -- Vt(bnd_out) = +Vt(bnd_in)
        if (msign_t.eq.+1) then 
          call bdsppy(mbnd,margint,vxm,ix,jx)
          call bdspny(mbnd,margint,vxdym,ix,jx)
          call bdsppy(mbnd,margint,vxdxm,ix,jx)
c       -- Vt(bnd_out) = -Vt(bnd_in)
        else
          call bdspny(mbnd,margint,vxm,ix,jx)
          call bdsppy(mbnd,margint,vxdym,ix,jx)
          call bdspny(mbnd,margint,vxdxm,ix,jx)
        endif

      endif
c     -- y-direction


      return
      end
