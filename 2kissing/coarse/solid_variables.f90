module solid_variables
  implicit none
  save
  
  real*8,parameter :: solid_scale = 0.5d0  !...scale the size of the structure
  integer,parameter :: nsurface = 4

  integer :: nn_solid  ,ne_solid
  integer :: nn_solid_1,ne_solid_1
  integer :: nsd_solid
  integer :: n_solid !...number of solids ( x times the solid, which is read in coortable
  integer :: nen_solid   !...number of nodes per element
  integer :: iquad_solid  !...quadratur type, see "quadxdxn.f"
  integer :: nquad_solid
  integer,parameter :: ndfpad_solid=5,nsdpad_solid=3,nenpad_solid=8,nquadpad_solid=8

  real* 8 xq_solid(nsdpad_solid,nquadpad_solid),wq_solid(nquadpad_solid)


  integer :: n_solid_ess_BC !...number of nodes with essential BC (displacement)

  real*8,dimension(:,:),allocatable :: solid_ess_BC
  real*8,allocatable :: shift(:,:)

!  real*8,dimension(:,:),allocatable,save :: solid_ess_BC
!  real*8,allocatable,save :: shift(:,:)

  integer :: mapping_solid(6,8,8)

contains




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solid_variables_BC_initialize
  implicit none
  
  integer :: error_id


  if (.not.allocated(solid_ess_BC)) then
     allocate(solid_ess_BC(1:nsd_solid,1:n_solid_ess_BC),stat=error_id)
  endif


  write(*,*) "memory for solid variables BC succesful allocated"

end subroutine solid_variables_BC_initialize




end module solid_variables
