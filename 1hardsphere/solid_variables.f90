module solid_variables
  implicit none
  save
  
  real*8,parameter :: solid_scale = 2.0d0  !...scale the size of the structure
  integer,parameter :: nsurface = 4

  integer :: nn_solid  ,ne_solid
  integer :: nn_solid_1,ne_solid_1
  integer :: nsd_solid
  integer :: n_solid !...number of solids ( x times the solid, which is read in coortable
  integer :: nis   !...number of nodes per element

  integer,dimension(:,:),allocatable :: solid_fem_con   !...connectivity for solid FEM mesh
  integer,dimension(:,:),allocatable :: solid_surface   !...surface element faces

  real*8,dimension(:,:),allocatable :: solid_force_FSI   !...fluid structure interaction force
  real*8,dimension(:,:),allocatable :: solid_coor_init   !...node position initial
  real*8,dimension(:,:),allocatable :: solid_coor_curr   !...node position current
  real*8,dimension(:,:),allocatable :: solid_vel         !...velocity
  real*8,dimension(:,:),allocatable :: solid_prevel      !...velocity - previous timestep
  real*8,dimension(:,:),allocatable :: solid_accel       !...acceleration


  real*8,dimension(:)  ,allocatable :: pave(:)  !...averaged solid pressure (from mixed formulation -> ???)

  real*8,dimension(:,:),allocatable :: solid_stress  !...solid stress (Voigt notation)
  real*8,dimension(:,:),allocatable :: solid_strain  !...solid strain (Voigt notation)

  integer :: n_solid_ess_BC !...number of nodes with essential BC (displacement)

  real*8,dimension(:,:),allocatable :: solid_ess_BC

  integer :: mapping_solid(6,8,8)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solid_variables_initialize
  use error_memory
  implicit none

  if (allocated(solid_fem_con))   deallocate(solid_fem_con)
  if (allocated(solid_surface))   deallocate(solid_surface)
  if (allocated(solid_force_FSI)) deallocate(solid_force_FSI)
  if (allocated(solid_coor_init)) deallocate(solid_coor_init)
  if (allocated(solid_coor_curr)) deallocate(solid_coor_curr)
  if (allocated(solid_vel))       deallocate(solid_vel)
  if (allocated(solid_prevel))    deallocate(solid_prevel)
  if (allocated(solid_accel))     deallocate(solid_accel)


  allocate(solid_fem_con(1:ne_solid,1:nis)        ,stat=error_id)
  allocate(solid_surface(1:ne_solid,1:nsurface)   ,stat=error_id)
  allocate(solid_force_FSI(1:nsd_solid,1:nn_solid),stat=error_id)
  allocate(solid_coor_init(1:nsd_solid,1:nn_solid),stat=error_id)
  allocate(solid_coor_curr(1:nsd_solid,1:nn_solid),stat=error_id)
  allocate(solid_vel(      1:nsd_solid,1:nn_solid),stat=error_id)
  allocate(solid_prevel(   1:nsd_solid,1:nn_solid),stat=error_id)
  allocate(solid_accel(    1:nsd_solid,1:nn_solid),stat=error_id)


  if (.not.allocated(solid_stress)) then
     allocate(solid_stress(6,nn_solid),stat=error_id)
  endif
    if (.not.allocated(solid_strain)) then
     allocate(solid_strain(6,nn_solid),stat=error_id)
  endif
    if (.not.allocated(pave)) then
     allocate(pave(nn_solid),stat=error_id)
  endif


  write(*,*) "memory for solid variables succesful allocated"

end subroutine solid_variables_initialize



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solid_variables_BC_initialize
  use error_memory
  implicit none



  if (.not.allocated(solid_ess_BC)) then
     allocate(solid_ess_BC(1:nsd_solid,1:n_solid_ess_BC),stat=error_id)
  endif


  write(*,*) "memory for solid variables BC succesful allocated"

end subroutine solid_variables_BC_initialize




end module solid_variables
