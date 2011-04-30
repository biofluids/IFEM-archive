

module  denmesh_variables
  implicit none
  save

  integer :: nn_den   !# of dense  mesh points
  integer :: ne_den   !# of dense  mesh elements
  integer :: nen_den
  integer :: nbc_den  !# of boundary nodes for dense mesh
end module denmesh_variables
