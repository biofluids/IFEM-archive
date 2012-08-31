subroutine readpartele_solid(partele)
  use solid_variables, only: ne_solid, ne_solid_1, n_solid
  implicit none

  integer :: partele(ne_solid)
  integer :: file,i
  integer j

  file=21
  open(file, FILE="partele_solid.in", STATUS="old")
  do i=1,ne_solid_1
     read(file,*) partele(i)
  enddo

  do j=2,n_solid
  	do i=1,ne_solid_1
	partele((j-1)*ne_solid_1 + i) = partele(i)
	end do
  end do

  close(file)

  return
end subroutine readpartele_solid


subroutine readpartnode_solid(partnode)
  use solid_variables, only: nn_solid, nn_solid_1, n_solid
  implicit none

  integer :: partnode(nn_solid)
  integer :: file,i
  integer j
  file=21
  open(file, FILE="partnode_solid.in", STATUS="old")
  do i=1,nn_solid_1
     read(file,*) partnode(i)
  enddo

  do j=2,n_solid
  	do i=1,nn_solid_1
	partnode((j-1)*nn_solid_1 + i) = partnode(i)
	end do
  end do

  close(file)

  return
end subroutine readpartnode_solid
