subroutine r_sinit(solid_coor_init,solid_coor_curr)
  use solid_variables, only: nsd_solid,nn_solid
  use r_common
  implicit none

  real*8,dimension(1:nsd_solid,1:nn_solid) :: solid_coor_init   !...node position initial
  real*8,dimension(1:nsd_solid,1:nn_solid) :: solid_coor_curr   !...node position current
!cccccccccccccccccccccccccc
!     read initial velocity
!cccccccccccccccccccccccccc
  if (ninit .eq. 1) then
     call r_sreadinit(solid_coor_init,solid_coor_curr)
  endif
  return
end	subroutine r_sinit
