c     ----------------------------------------------------------------
      subroutine slvsml(u,rhs)
      double precision rhs(-1:2,-1:2),u(-1:2,-1:2)
      double precision hh
      hh=0.5d0
c     --- begin: periodic boundary condition
      u(0,0)=0.0d0
      u(1,1)=u(0,0)+hh*hh*(rhs(0,0)-rhs(1,1))/4.0d0
      u(1,0)=u(0,0)+hh*hh*(rhs(0,0)-rhs(1,1)-2.0d0*rhs(1,0))/8.0d0
      u(0,1)=u(0,0)+hh*hh*(rhs(0,1)-rhs(1,1)-2.0d0*rhs(1,0))/8.0d0
      u(2,0)=u(0,0)
      u(0,2)=u(0,0)
      u(2,2)=u(0,0)
      u(1,2)=u(1,0)
      u(2,1)=u(0,1)
c     --- end: periodic boundary condition
      return
      end
