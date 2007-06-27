subroutine r_sreadinit(solid_coor_init,solid_coor_curr)
  use r_common
  use solid_variables
  implicit none

  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_init   !...node position initial
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_curr   !...node position current
	integer:: i

!USER'S DEPENDENT
do i=1,nn_solid
	!solid_coor_curr(2,i) = (0.375*solid_coor_init(1,i)+0.925)*solid_coor_init(2,i)
	!solid_coor_curr(2,i) = (-0.3125*(solid_coor_init(1,i))**2+ &
	!		0.625*solid_coor_init(1,i)+0.9875)*solid_coor_init(2,i)
	!solid_coor_curr(2,i) = (-0.625*(solid_coor_init(1,i))**2+ &
	!		1.25*solid_coor_init(1,i)+0.675)*solid_coor_init(2,i)
enddo


  return
end subroutine r_sreadinit