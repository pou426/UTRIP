c======================================================================|
      function rangen(mseed)
c======================================================================|
c     random number generation
c     based on 
c        Park and Miller (1988, Comm. ACM, 31, 1192)
c        Schrage (1979, ACM Trans. Math. Soft., 5, 132)
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      ma=16807
      mm=2147483647
      mq=127773
      mr=2836
      rnml=4.656612875245796924105750827168d-10

      k=mseed/mq
      mseed=ma*(mseed-k*mq)-mr*k
      if (mseed.lt.0) mseed=mseed+mm
      rangen=rnml*mseed

      return
      end
