subroutine r_sinit
  use r_common
  implicit none
!cccccccccccccccccccccccccc
!     read initial velocity
!cccccccccccccccccccccccccc
  if (ninit .eq. 1) then
     call r_sreadinit
  endif
  return
end	subroutine r_sinit
