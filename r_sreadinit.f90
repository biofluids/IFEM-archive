subroutine r_sreadinit(solid_coor_init,solid_coor_curr)
  use r_common
  use solid_variables
  implicit none

  real*8,dimension(1:nsd_solid,1:nn_solid) :: solid_coor_init   !...node position initial
  real*8,dimension(1:nsd_solid,1:nn_solid) :: solid_coor_curr   !...node position current

  integer :: ndum,ndumtest
  real*8 :: x1,y1,z1

!...assign nonlinear initial conditions
 1000 read(1,*) ndum,x1,y1,z1
  ndumtest=ndum-1
  if (ndumtest .ge. 0) then
     if (initdir .eq. 1) then
        solid_coor_curr(1,ndum) = solid_coor_init(1,ndum) + x1
		solid_coor_curr(2,ndum) = solid_coor_init(2,ndum) + y1
		solid_coor_curr(3,ndum) = solid_coor_init(3,ndum) + z1

        !xindis(1,ndum)=x1
        !xindis(2,ndum)=y1
        !xindis(3,ndum)=z1
     endif
     goto 1000
  endif

  return
end	subroutine r_sreadinit