c     
c     caculation of fp
c     
      subroutine r_scalfp(fp,ocpp,i)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      fp=-ocpp*(bpre-cpre)*hp(i)
      return
      end


