      subroutine r_sinit
      implicit real*8 (a-h,o-z)
      include 'r_common'
ccccccccccccccccccccccccccc
c     read initial velocity
ccccccccccccccccccccccccccc
      if (ninit .eq. 1) then
         call r_sreadinit
      endif
      return
      end	