

module  centermesh_variables
  implicit none
  save


  integer :: flag_center !0 for matlab 1 for mesh info
  integer :: nn_center   !# of center mesh points
  integer :: ne_center   !# of center mesh elements
  integer :: nen_center
end module centermesh_variables
