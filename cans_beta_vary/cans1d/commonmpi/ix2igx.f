      subroutine ix2igx(ix,margin
     &  ,igx,ig0,ig1
     &  ,mpex,ipex,mpe,ipe,mpex0,merr)

      if (mpe.eq.0) then
        write(*,*) 'ipex,mpex,mpex0= '
     &        , ipex, mpex, mpex0
      endif
      if (ipex.eq.0) then
         write(*,*) 'ERROR in mpe '
         merr=1
         return
      endif
      if (mpex.ne.mpex0) then
         write(*,*) 'Mismatch in # of PEs !'
         write(*,*) 'ipex,mpex,mpex0= '
     &        , ipex, mpex, mpex0
         merr=1
         return
      endif

      igx=ix*ipex-2*margin*(ipex-1)
      if (mpe.eq.0) then
        write(6,*) 'igx= ',igx
      endif

      ipe=mpe

      ig0=ipe*(ix-2*margin)+1+margin
      ig1=ig0+(ix-2*margin)-1

      return
      end
