module solid_fem_BC
  implicit none
  save

  integer,private :: nBC_ess_type_id
  integer,dimension(:)  ,allocatable,private :: BC_ess_type_id
  real*8 ,dimension(:,:),allocatable,private :: BC_ess_type

  integer,private :: nn_solid_BC_ess
  integer,dimension(:)  ,allocatable,private :: solid_BC_ess
  integer,dimension(:,:),allocatable,private :: solid_BC_ess_value

  integer,parameter,private :: typefile_unit = 20
  character(len=27),parameter,private :: typefile = "input_solid_BC_ess_types.in"


 !...private subroutines
  private :: read_BC_ess_types


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solid_fem_BC_apply_essential
!
! Axel Gerstenberger, NWU, Mai 2003
!
! This subroutine applies "penalty" forces on nodes with Dirichlet (essential) Boundary Conditions
! 
  use solid_variables, only: solid_force_FSI,solid_coor_init,solid_coor_curr
  implicit none

  integer :: innBC

  real*8,dimension(1:3) :: xx,yy,dist,force_BC

  write(*,*) " apply essential BC "

  do innBC = 1,nn_solid_BC_ess
     xx(1:3) = solid_coor_init(1:3,solid_BC_ess(innBC))
	 yy(1:3) = solid_coor_curr(1:3,solid_BC_ess(innBC))

     dist(1:3)         = yy(1:3) - xx(1:3)

	 force_BC(1:3) = dist(1:3) * solid_BC_ess_value(1:3,innBC)

	 solid_force_FSI(1:3,solid_BC_ess(innBC)) = solid_force_FSI(1:3,solid_BC_ess(innBC)) - force_BC(1:3)
  enddo

end subroutine solid_fem_BC_apply_essential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solid_fem_BC_read_essential
  use solid_variables
  implicit none

  integer :: test_node_BC_id(1:nn_solid)
  real*8 :: test_node_BC_value(1:3,1:nn_solid)

  integer,parameter :: ifileunit=20

  integer :: numgb,i,j,inn,innBC,test_node_BC_type,nodes_per_BC_type,n_test_node
  integer :: counter_solid_BC_node,BC_ess_type_pos

  call read_BC_ess_types

  test_node_BC_id(:) = 0
  test_node_BC_value(1:3,1:nn_solid) = 0.0d0

  open(unit = ifileunit,file = "input_solid_BC.in",status="old",action="read")

  read(ifileunit,*) numgb
  write(*,*) "number of different BC types defined = ",numgb

  do i=1,numgb  
     read(ifileunit,*) test_node_BC_type   !...ndirgb (old)
	 write(*,*) i,". BC:  ",test_node_BC_type
	 do j = 1,nBC_ess_type_id  
        if (BC_ess_type_id(j) == test_node_BC_type)  BC_ess_type_pos = j !...find vector index for BC_type
		!write(*,*) BC_ess_type_pos
	 enddo
     
	 read(ifileunit,*) nodes_per_BC_type   !...numdir (old)
	 write(*,*)   "  ",nodes_per_BC_type

	 do innBC = 1,nodes_per_BC_type
	    read(ifileunit,*) n_test_node
		write(*,*) "node:",n_test_node
        test_node_BC_id(n_test_node) = 1
		test_node_BC_value(1:3,n_test_node) = test_node_BC_value(1:3,n_test_node) + BC_ess_type(1:3,BC_ess_type_pos)
	 enddo
  enddo 

  close(ifileunit)


  nn_solid_BC_ess = sum(test_node_BC_id(:))  !...number of nodes that have a BC applied

  allocate(solid_BC_ess(1:nn_solid_BC_ess))            !...contains node number
  allocate(solid_BC_ess_value(1:3,1:nn_solid_BC_ess))  !...contains three penalty parameter for force calculation
  
  counter_solid_BC_node = 0
  do inn = 1,nn_solid
     if (test_node_BC_id(inn) == 1) then  !...if BC applied
	    counter_solid_BC_node = counter_solid_BC_node + 1
	    solid_BC_ess(counter_solid_BC_node) = inn
        solid_BC_ess_value(1:3,counter_solid_BC_node) = test_node_BC_value(1:3,inn)
	 endif
  enddo
  

  write(*,*) "BC for ",nn_solid_BC_ess," nodes sucessfully read..."

end subroutine solid_fem_BC_read_essential


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_BC_ess_types
  implicit none

  integer :: i

  open(unit=typefile_unit,file=typefile,status="old",action="read")

  read(typefile_unit,*) nBC_ess_type_id
  write(*,*) "number of boundary types defined",nBC_ess_type_id

  allocate(BC_ess_type_id( 1:nBC_ess_type_id))
  allocate(BC_ess_type(1:3,1:nBC_ess_type_id))

  do i = 1,nBC_ess_type_id
     read(typefile_unit,*) BC_ess_type_id(i) , BC_ess_type(1:3,i)
  enddo

  !write(*,*) BC_ess_type

  close(typefile_unit)

end subroutine read_BC_ess_types



end module solid_fem_BC