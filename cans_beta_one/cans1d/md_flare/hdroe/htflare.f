c***********************************************************************
      subroutine htflare(pr,dt,gm,t,x,ix)
c***********************************************************************
c
c  heating function of flare
c
c***********************************************************************
      implicit double precision (a-h,o-z)

      dimension x(ix)
      dimension pr(ix)
      dimension flr(ix)
c----------------------------------------------------------------------|
c  flare condition
c----------------------------------------------------------------------|
c  dzflr : distance from loop top to flare site
      pi = acos(-1.0d0)

      dzflr=0.d0
      wflr=30.d0
      zlim=20.d0
      zflr=x(ix-1)-dzflr
c----------------------------------------------------------------------|

      p2r=sqrt(2*pi)

      do i=1,ix
         flr(i)=exp(-0.5d0*((x(i)-zflr)/wflr)**2)/p2r
     &          *0.5d0*(tanh((x(i)-zlim)/3.d0)+1.)
      enddo
         flr(ix)=flr(ix-2)
c======================================================================@
c  flare condition
c======================================================================@
c   energy input to the half loop 
c     en_flare = ein*(tinend-tinstt)
c                *rkb*tenml*denml*rlnml [erg/cm^3 * cm]
c     en_flare = 1e12 [erg/cm^2] * (ein/0.02)*(tinend/12)
c                *(tenml/1.e4K)*(denml/1.e17cm^-3)*(rlnml/2.e7cm)
c
      ein=0.0005d0
      tinstt=0.d0
      tinend=12.d0
      wstt=0.1d0
      wend=0.1d0
c----------------------------------------------------------------------c
c  flare heating
c----------------------------------------------------------------------c

      q0=(  tanh((t-tinstt)/wstt) + 1.)*0.5d0
      q1=( -tanh((t-tinend)/wend) + 1.)*0.5d0
      q=q0*q1

      do i=1,ix
        pr(i)=pr(i)+dt*(gm-1)/gm*ein*q*flr(i)
      enddo


      return
      end
