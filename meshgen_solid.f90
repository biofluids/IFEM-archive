module meshgen_solid
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readx_solid(xyz,nn,nsd)
  implicit none

  integer,intent(in) :: nn,nsd
  real(8) :: xyz(nn,nsd)
  integer :: idummy,inn,file

  file=23
  open(file, FILE="mxyz_solid.in", STATUS="old",action="read")

  do inn=1,nn
     read(file,*) idummy,xyz(inn,1:nsd)
  enddo

  close(file)

  return
end subroutine readx_solid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readien_solid(solid_con,ne,nen)
  implicit none

  integer,intent(in) :: ne,nen
  integer :: solid_con(ne,nen)

  integer :: file,ine,idummy

  file=21
  open(file, FILE="mien_solid.in", STATUS="old",action="read")

  do ine=1,ne
     read(file,*) idummy,solid_con(ine,1:nen),idummy
  enddo
  close(file)

  return
end subroutine readien_solid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!******************
! Lucy commented this out, because we are not using it.
!******************
!subroutine readrng_solid(rngface,ne,neface)
!  implicit none

!  integer,intent(in) :: ne,neface
!  integer :: rngface(ne,neface)

!  integer :: file,i,ieface,iec

!  file=21
!  open(file, FILE="mrng_solid.in", STATUS="old",action="read")
!  do i=1,ne
!     read(file,*) rngface(i,:)
!  enddo

!  do ieface=1,neface
!     do iec=1,ne
!        if(rngface(iec,ieface).lt.0) rngface(iec,ieface) = 0
!     enddo
!  enddo

!  close(file)

!  return
!end subroutine readrng_solid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_solid_ale_boundary(solid_ale_boundary)
  use solid_variables, only: nn_solid
  implicit none

  integer,intent(out) :: solid_ale_boundary(nn_solid)

  integer :: file,innBC,n_test_node,n_ale_boundary

  file=21
  open(file, FILE="input_solid_ale_boundary.in", STATUS="old",action="read")

  solid_ale_boundary = 0

  read(file,*) n_ale_boundary
  do innBC = 1,n_ale_boundary
     read(file,*) n_test_node
     write(*,*) innBC,". node:",n_test_node
     solid_ale_boundary(n_test_node) = 1
  enddo

  close(file)

  return
end subroutine read_solid_ale_boundary


end module meshgen_solid