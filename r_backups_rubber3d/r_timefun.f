c     time function
      subroutine r_timefun
      implicit real*8 (a-h,o-z)
      include 'r_common'
      pi=datan(1.0d0)*4.0d0
ccccccccccccccccccccccccccccccccc
c     timefunction type 3
ccccccccccccccccccccccccccccccccc
      if (ntfun .eq. 3) then
c         if (iti .eq. 4) then
c            tfun(3)=1.0d0
c         else
c            tfun(3)=0.0d0
c         endif
c         tfun(3)=1.0d0*dsin(2.0d0*pi*iti/50)
         tfun(3)=1.0d0*dsin(pi*(iti-1)/2.0d0)
c         tfun(3)=1.0d0
      endif
ccccccccccccccccccccccccccccccccc      
c     timefunction type 2
ccccccccccccccccccccccccccccccccc
      if (ntfun .eq. 2) then
         tfun(2)=0.1d0
      endif
ccccccccccccccccccccccccccccccccc      
c     timefunction type 1
ccccccccccccccccccccccccccccccccc
      if (ntfun .eq. 1) then
         tfun(1)=1.0d0*iti/nts_solid
      endif
ccccccccccccccccccccccccccccccccc      
c     timefunction type 4
ccccccccccccccccccccccccccccccccc
      if (ntfun .eq. 4) then
         tfun(4)=1.0d0
      endif
ccccccccccccccccccccccccccccccccc      
c     timefunction type 5
ccccccccccccccccccccccccccccccccc
      if (ntfun .eq. 5) then
         if (iti .ge. 90) then
            tfun(5)=1.0d0*dsin(2.0d0*pi*(iti-90)/70.0d0)
         else
            tfun(5)=0.0d0
         endif
      endif
c
      return
      end
