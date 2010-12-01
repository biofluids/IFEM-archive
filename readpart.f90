subroutine readpartele(partele)
  use fluid_variables, only: ne
  implicit none

  integer :: partele(ne)
  integer :: file,i

  file=21
  open(file, FILE="partele.in", STATUS="old")
  do i=1,ne
     read(file,*) partele(i)
  enddo
  close(file)

  return
end subroutine readpartele


subroutine readpartnode(partnode)
  use fluid_variables, only: nn
  implicit none

  integer :: partnode(nn)
  integer :: file,i

  file=21
  open(file, FILE="partnode.in", STATUS="old")
  do i=1,nn
     read(file,*) partnode(i)
  enddo
  close(file)

  return
end subroutine readpartnode
