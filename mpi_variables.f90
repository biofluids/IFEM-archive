module mpi_variables
implicit none
save
integer ierror
integer ncpus    
! ncpus = count of cores used, say 1024 as partitioned and defined when actually when running cases
! commented by Jubiao Yang on 03/14/2013
integer myid
end module mpi_variables
