      subroutine ix2igx(ix,jx,kx,margin
     &  ,igx,jgx,kgx,ig0,ig1,jg0,jg1,kg0,kg1
     &  ,mpex,ipex,jpex,kpex,mpe,ipe,jpe,kpe,mpex0,merr)

      if (mpe.eq.0) then
        write(*,*) 'ipex,jpex,kpex,mpex,mpex0= '
     &        , ipex, jpex, kpex, mpex, mpex0
      endif
      if ( (ipex.eq.0).or.(jpex.eq.0).or.(kpex.eq.0) )then
         write(*,*) 'ERROR in mpe '
         merr=1
         return
      endif
      if (mpex.ne.mpex0) then
         write(*,*) 'Mismatch in # of PEs !'
         write(*,*) 'ipex,jpex,kpex,mpex,mpex0= '
     &        , ipex, jpex, kpex, mpex, mpex0
         merr=1
         return
      endif

      igx=ix*ipex-2*margin*(ipex-1)
      jgx=jx*jpex-2*margin*(jpex-1)
      kgx=kx*kpex-2*margin*(kpex-1)
      if (mpe.eq.0) then
        write(6,*) 'igx,jgx,kgx= ',igx,jgx,kgx
      endif

      kpe=mpe/(ipex*jpex)
      jpe=(mpe-ipex*jpex*kpe)/ipex
      ipe=mpe-ipex*jpex*kpe-ipex*jpe

      ig0=ipe*(ix-2*margin)+1+margin
      ig1=ig0+(ix-2*margin)-1
      jg0=jpe*(jx-2*margin)+1+margin
      jg1=jg0+(jx-2*margin)-1
      kg0=kpe*(kx-2*margin)+1+margin
      kg1=kg0+(kx-2*margin)-1

      return
      end
