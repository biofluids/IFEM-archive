module solid_variables
  implicit none
  save
  
  public

  real*8,parameter :: solid_scale = 1.0d0  !...scale the size of the structure

  integer :: nn_solid  ,ne_solid
  integer :: nn_solid_1,ne_solid_1
  integer :: nsd_solid
  integer :: n_solid !...number of solids ( x times the solid, which is read in coortable

  real*8,dimension(:,:),allocatable :: solid_force_FSI   !...fluid structure interaction force
  real*8,dimension(:,:),allocatable :: solid_coor_init   !...node position initial
  real*8,dimension(:,:),allocatable :: solid_coor_curr   !...node position current
  real*8,dimension(:,:),allocatable :: solid_vel         !...velocity
  real*8,dimension(:,:),allocatable :: solid_prevel      !...velocity - previous timestep
  real*8,dimension(:,:),allocatable :: solid_accel       !...acceleration


  integer :: n_solid_ess_BC !...number of nodes with essential BC (displacement)

  real*8,dimension(:,:),allocatable :: solid_ess_BC

contains

subroutine solid_variables_initialize
  implicit none
  integer :: err

	write(*,*) 'fine here'
  if (.not.allocated(solid_force_FSI)) then
     allocate(solid_force_FSI(1:nsd_solid,1:nn_solid),stat=err)
  endif
  if (.not.allocated(solid_coor_init)) then
     allocate(solid_coor_init(1:nsd_solid,1:nn_solid),stat=err)
  endif
  if (.not.allocated(solid_coor_curr)) then
     allocate(solid_coor_curr(1:nsd_solid,1:nn_solid),stat=err)
  endif
  if (.not.allocated(solid_vel)) then
     allocate(solid_vel(      1:nsd_solid,1:nn_solid),stat=err)
  endif
  if (.not.allocated(solid_prevel)) then
     allocate(solid_prevel(   1:nsd_solid,1:nn_solid),stat=err)
  endif
  if (.not.allocated(solid_accel)) then
     allocate(solid_accel(    1:nsd_solid,1:nn_solid),stat=err)
  endif


  write(*,*) "memory for solid variables succesful allocated"

end subroutine solid_variables_initialize




subroutine solid_variables_BC_initialize
  implicit none
  integer :: err


  if (.not.allocated(solid_ess_BC)) then
     allocate(solid_ess_BC(1:nsd_solid,1:n_solid_ess_BC),stat=err)
  endif


  write(*,*) "memory for solid variables BC succesful allocated"

end subroutine solid_variables_BC_initialize




end module solid_variables
