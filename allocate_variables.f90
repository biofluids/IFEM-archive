!=============================================
!used for allocatable variables
!=============================================

module allocate_variables
  implicit none
  save

  integer,dimension(:),allocatable :: den_domain
  integer,dimension(:),allocatable :: center_domain
  integer :: ne_den_domain
  integer :: nn_center_domain

  integer,dimension(:),allocatable :: inter_ele
  integer,dimension(:),allocatable :: inter_ele_den
  integer,dimension(:),allocatable :: regen_ele
  integer :: ne_inter,ne_inter_den,ne_regen_ele
  integer,dimension(:),allocatable :: regen_ele_loc
  integer :: ne_regen_ele_loc
end module allocate_variables
