c     ----------------------------------------------------------------
      subroutine melloc(len,qq,ng,ngrid,z,memlen,mem)
      integer len,qq(ng),ng,ngrid,memlen,mem
      double precision z(memlen)
      if (mem+len+1.gt.memlen) pause 'insufficient memory in melloc'
      z(mem+1)=len
      qq(ngrid)=mem+2
      mem=mem+len+1
      return
      end
