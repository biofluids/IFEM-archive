c     
c     caculation of kpp
c     
      subroutine r_scalkpp(fkpp,ocpp,k,m)
      implicit real*8 (a-h,o-z)
      include 'r_common'
      fkpp=ocpp*hp(k)*hp(m)
      return
      end
