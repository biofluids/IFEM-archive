subroutine r_timefun
  use global_constants
  use r_common
  implicit none

!ccccccccccccccccccccccccccccccccc
!     timefunction type 3
!ccccccccccccccccccccccccccccccccc
  if (ntfun == 3) then
!         if (iti == 4) then
!            tfun(3)=1.0d0
!         else
!            tfun(3)=0.0d0
!         endif
!         tfun(3)=1.0d0*dsin(2.0d0*pi*iti/50)
!         tfun(3)=1.0d0*dsin(pi*(iti-1)/2.0d0)
!         tfun(3)=1.0d0
  endif
!ccccccccccccccccccccccccccccccccc      
!     timefunction type 2
!ccccccccccccccccccccccccccccccccc
  if (ntfun == 2) then
     tfun(2)=0.1d0
  endif
!ccccccccccccccccccccccccccccccccc      
!     timefunction type 1
!ccccccccccccccccccccccccccccccccc
!      if (ntfun == 1) then
!         tfun(1)=1.0d0*iti/nts_solid
!      endif
!ccccccccccccccccccccccccccccccccc      
!     timefunction type 4
!ccccccccccccccccccccccccccccccccc
  if (ntfun == 4) then
     tfun(4)=1.0d0
  endif
!ccccccccccccccccccccccccccccccccc      
!     timefunction type 5
!ccccccccccccccccccccccccccccccccc
!      if (ntfun == 5) then
!         if (iti >= 90) then
!            tfun(5)=1.0d0*dsin(2.0d0*pi*(iti-90)/70.0d0)
!         else
!            tfun(5)=0.0d0
!         endif
!      endif
!
  return
end subroutine r_timefun
