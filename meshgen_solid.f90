module meshgen_solid
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readx_solid(xyz,nn,nsd)
  implicit none

  integer,intent(in) :: nn,nsd
  real(8) :: xyz(nn,nsd)
  integer :: inn,file

  file=23
  open(file, FILE="mxyz_solid.in", STATUS="old",action="read")

  do inn=1,nn
     read(file,*) xyz(inn,1:nsd)
  enddo

  close(file)

  return
end subroutine readx_solid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readien_solid(solid_con,ne,nen)
  implicit none

  integer,intent(in) :: ne,nen
  integer :: solid_con(ne,nen)

  integer :: file,ine

  file=21
  open(file, FILE="mien_solid.in", STATUS="old",action="read")

  do ine=1,ne
     read(file,*) solid_con(ine,1:nen)
  enddo
  close(file)

  return
end subroutine readien_solid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_sfcon(sfcon,node_sfcon)
  implicit none

  integer,intent(in) :: node_sfcon
  integer :: sfcon(node_sfcon)
  integer :: file,i

  file=23
  open(file, FILE="sfcon.in", STATUS="old",action="read")

  do i=1,node_sfcon
     read(file,*) sfcon(i)
  enddo

  close(file)

  return
end subroutine read_sfcon


end module meshgen_solid
