      subroutine ix2igx(ix,jx,margin
     &  ,igx,jgx,ig0,ig1,jg0,jg1
     &  ,mpex,ipex,jpex,mpe,ipe,jpe,mpex0,merr)

      if (mpe.eq.0) then
        write(*,*) 'ipex,jpex,mpex,mpex0= '
     &        , ipex, jpex, mpex, mpex0
      endif
      if ( (ipex.eq.0).or.(jpex.eq.0) )then
         write(*,*) 'ERROR in mpe '
         merr=1
         return
      endif
      if (mpex.ne.mpex0) then
         write(*,*) 'Mismatch in # of PEs !'
         write(*,*) 'ipex,jpex,mpex,mpex0= '
     &        , ipex, jpex, mpex, mpex0
         merr=1
         return
      endif

      igx=ix*ipex-2*margin*(ipex-1)
      jgx=jx*jpex-2*margin*(jpex-1)
      if (mpe.eq.0) then
        write(6,*) 'igx,jgx= ',igx,jgx
      endif

      jpe=mpe/ipex
      ipe=mpe-ipex*jpe

      ig0=ipe*(ix-2*margin)+1+margin
      ig1=ig0+(ix-2*margin)-1
      jg0=jpe*(jx-2*margin)+1+margin
      jg1=jg0+(jx-2*margin)-1

      return
      end
