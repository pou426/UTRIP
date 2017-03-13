c======================================================================|
      subroutine cipbdbsym(bxm,bym
     &                    ,margin,ix,jx,mdir,mbnd,msign_n,msign_t)
c======================================================================|
c----------------------------------------------------------------------|
      implicit none
      integer margin,ix,jx
      integer mdir,mbnd,msign_n,msign_t
      integer marginn,margint
      double precision bxm
      double precision bym
      dimension bxm(ix,jx),bym(ix,jx)
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

c       -- Bn(bnd_out) = +Bn(bnd_in)
        if (msign_n.eq.+1) then 
          call bdsmpx(mbnd,marginn,bxm,ix,jx)
c       -- Bn(bnd_out) = -Bn(bnd_in)
        else
          call bdsmnx(mbnd,marginn,bxm,ix,jx)
        endif

c       -- Bt(bnd_out) = +Bt(bnd_in)
        if (msign_t.eq.+1) then 
          call bdsppx(mbnd,margint,bym,ix,jx)
c       -- Bt(bnd_out) = -Bt(bnd_in)
        else
          call bdspnx(mbnd,margint,bym,ix,jx)
        endif

      endif
c     -- x-direction

c     -- y-direction
      if (mdir.eq.1) then 

c       -- Bn(bnd_out) = +Bn(bnd_in)
        if (msign_n.eq.+1) then 
          call bdsmpy(mbnd,marginn,bym,ix,jx)
c       -- Bn(bnd_out) = -Bn(bnd_in)
        else
          call bdsmny(mbnd,marginn,bym,ix,jx)
        endif

c       -- Bt(bnd_out) = +Bt(bnd_in)
        if (msign_t.eq.+1) then 
          call bdsppy(mbnd,margint,bxm,ix,jx)
c       -- Bt(bnd_out) = -Bt(bnd_in)
        else
          call bdspny(mbnd,margint,bxm,ix,jx)
        endif

      endif
c     -- y-direction


      return
      end
