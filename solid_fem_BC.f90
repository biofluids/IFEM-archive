module solid_fem_BC
!  use solid_variables, only:nn_solid,nsd_solid
! integer :: solid_BC_ess(nsd_solid)
!  integer :: solid_BC_ess_value(nsd_solid,nn_solid)

!  integer :: nBC_ess_type_id
!  integer :: BC_ess_type_id(20)
!  real(8) :: BC_ess_type(nsd_solid,20)

  implicit none
  save

!  integer,private :: nBC_ess_type_id
!  use solid_variables, only:nn_solid,nsd_solid
!  integer :: BC_ess_type_id(20)
!  real(8) :: BC_ess_type(nsd_solid,20)
!  integer :: nn_solid_BC_ess
!  integer:: nn_solid,nsd_solid
!  use solid_variables, only:nn_solid,nsd_solid
  integer:: solid_BC_ess(3)
  integer:: solid_BC_ess_value(3,1:1000000)  !this is a temporary solution, needs to figure out how to use other modules variables here

!  integer,parameter,private :: typefile_unit = 20
!  character(len=27),parameter,private :: typefile = "input_solid_BC_ess_types.in"


 !...private subroutines
!  private :: read_BC_ess_types


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine read_BC_ess_types
!  use solid_variables, only: nsd_solid
  
!  implicit none

!  integer :: i,error_id
!  integer,private :: nBC_ess_type_id
!  integer,private :: BC_ess_type_id(20)
!  real(8),private :: BC_ess_type(nsd_solid,20)
!  integer,parameter,private :: typefile_unit = 20
!  character(len=27),parameter,private :: typefile = "input_solid_BC_ess_types.in"

!  open(unit=typefile_unit,file=typefile,status="old",action="read")

!  read(typefile_unit,*) nBC_ess_type_id
!  write(*,*) "number of boundary types defined",nBC_ess_type_id
!
!  if (allocated(BC_ess_type_id)) then
!     deallocate(BC_ess_type_id)
!  endif
!  if (allocated(BC_ess_type)) then
!     deallocate(BC_ess_type)
!  endif

!  allocate(BC_ess_type_id(nBC_ess_type_id),stat=error_id)
!  allocate(BC_ess_type(nsd_solid,nBC_ess_type_id),stat=error_id)

!  do i = 1,nBC_ess_type_id
!     read(typefile_unit,*) BC_ess_type_id(i) , BC_ess_type(1:nsd_solid,i)
!  enddo


 ! close(typefile_unit)

!end subroutine read_BC_ess_types

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solid_fem_BC_read_essential(n_bc,bc_type)
  use solid_variables, only: nsd_solid,nn_solid
  implicit none

  integer,parameter :: ifileunit=20
  integer :: i,inn,innBC,bc_type_ID
  integer,dimension(1:6,1:nn_solid) :: bc_type   !...list of boundary condition types per node
  integer,dimension(1:nn_solid) :: n_bc   !...number of boundary conditions per node
!  integer :: solid_BC_ess(nsd_solid)
!  integer :: solid_BC_ess_value(nsd_solid,nn_solid)

  integer :: nBC_ess_type_id
  integer :: BC_ess_type_id(20)
  real(8) :: BC_ess_type(nsd_solid,20)
  integer,parameter :: typefile_unit = 20
  character(len=27),parameter :: typefile = "input_solid_BC_ess_types.in"

  open(unit=typefile_unit,file=typefile,status="old",action="read")

  read(typefile_unit,*) nBC_ess_type_id
  write(*,*) "number of boundary types defined",nBC_ess_type_id

  do i = 1,nBC_ess_type_id
     read(typefile_unit,*) BC_ess_type_id(i) , BC_ess_type(1:nsd_solid,i)
  enddo


  close(typefile_unit)
 !call read_BC_ess_types
  write(*,*) 'okay'

  solid_BC_ess_value(:,:)=0.0
  do inn=1,nn_solid
    do innBC=1, n_bc(inn)
       bc_type_ID=bc_type(innBC,inn)
       solid_BC_ess_value(:, inn)=solid_BC_ess_value(:,inn)+BC_ess_type(:, bc_type_ID)
    enddo
  enddo
 write(*,*) 'BC_ess_type(:,6)=',BC_ess_type(:,6)
 write(*,*) 'solid_BC_ess_value(:,24)=',solid_BC_ess_value(:,24)
end subroutine solid_fem_BC_read_essential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solid_fem_BC_apply_essential(solid_force_FSI,solid_coor_init,solid_coor_curr)
!
! Axel Gerstenberger, NWU, Mai 2003
!
! This subroutine applies "penalty" forces on nodes with Dirichlet (essential) Boundary Conditions
! 
  use solid_variables, only: nn_solid,nsd_solid
  implicit none

  real(8),dimension(1:nsd_solid,1:nn_solid),intent(inout) :: solid_force_FSI   !...fluid structure interaction force
  real(8),dimension(1:nsd_solid,1:nn_solid),intent(in) :: solid_coor_init   !...node position initial
  real(8),dimension(1:nsd_solid,1:nn_solid),intent(in) :: solid_coor_curr   !...node position current
!  real(8),dimension(1:nsd_solid,1:nn_solid) ::solid_BC_ess_value            !...essential force values
  integer :: innBC

  real(8),dimension(1:nsd_solid) :: xx,yy,dist,force_BC

  write(*,*) " apply essential BC "
!  write(*,*) 'nn_solid_BC_ess=',nn_solid_BC_ess
!  write(*,*) 'solid_BC_ess_value=', solid_BC_ess_value(1:3,6)
!  do innBC = 1,nn_solid_BC_ess  !Lucy commented it out for new BC input format
   do innBC = 1, nn_solid
     xx(1:nsd_solid) = solid_coor_init(1:nsd_solid,innBC)
     yy(1:nsd_solid) = solid_coor_curr(1:nsd_solid,innBC)

     dist(1:nsd_solid)         = yy(1:nsd_solid) - xx(1:nsd_solid)
     force_BC(1:nsd_solid)=solid_BC_ess_value(1:nsd_solid,innBC)*0.5
     solid_force_FSI(1:nsd_solid,innBC) = solid_force_FSI(1:nsd_solid,innBC) - force_BC(1:nsd_solid)
  enddo

end subroutine solid_fem_BC_apply_essential



end module solid_fem_BC
