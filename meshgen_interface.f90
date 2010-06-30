module meshgen_interface
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readinter(xyz)
  use interface_variables, only: nn_inter
  use fluid_variables, only: nsd
  implicit none

  real(8) :: xyz(nsd,nn_inter)
  integer :: file,i

  file=21
  open(file,FILE='xyz_interface.in',status='old')

  do i=1,nn_inter
     read(file,*) xyz(1:nsd,i)
  end do

  close(file)
 
  return
end subroutine readinter

end module meshgen_interface
