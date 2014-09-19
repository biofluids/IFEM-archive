module mpi_variables
    implicit none
    save
    integer ierror
    integer ncpus
    integer myid
    
    integer, allocatable, dimension(:,:) :: sub_address
    ! sub_address(proc id, number of nodes share by this proc)
    integer, allocatable, dimension(:,:) :: sub_address_solid
    ! sub_address(proc id, number of nodes share by this proc) for solid mesh
    integer countrow
    ! length of sub_address
    integer countrow_solid
    ! length of sub_address_solid

end module mpi_variables