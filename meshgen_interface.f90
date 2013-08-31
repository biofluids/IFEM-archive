module meshgen_interface
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readinter(xyz)
  use interface_variables, only: nn_inter,maxmatrix
  use fluid_variables, only: nsd
  use mpi_variables
  implicit none

  real(8) :: xyz(nsd,maxmatrix)
  integer :: file,i

  file=21
  open(file,FILE='xyz_interface.in',status='old')
  if(nsd==2) then
    do i=1,nn_inter
       read(file,111) xyz(1,i),xyz(2,i)
    end do
111 format(f14.10,f14.10)
  else if( nsd==3) then
    do i=1,nn_inter
	read(file,112) xyz(1,i),xyz(2,i),xyz(3,i)
    end do
  end if
112 format(f14.10,f14.10,f14.10)
  close(file)


  return

end subroutine readinter

end module meshgen_interface
