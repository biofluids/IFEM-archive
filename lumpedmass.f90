module lumpedmass
  implicit none
  save

  real(8),allocatable :: mpqlumpmass(:)
  real(8),allocatable :: hattauij(:,:,:)  ! \hat{\tau_{ij}} and \R{\tau_{ij}}
                                   ! \hat{\tau_{ij}} = \R{\tau_{ij}} / M_A
  integer,allocatable :: flagdivision(:)
end module lumpedmass
