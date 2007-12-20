module meshgen_solid
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! added bc_type, n_bc in the subroutine
!-------------------------------------------------------
subroutine readx_solid(xyz,nn,nsd,n_bc,bc_type)
  implicit none

  integer,intent(in) :: nn,nsd
  real(8) :: xyz(nn,nsd)
  integer :: inn,file,n_bc(nn),bc_type(6,nn)

  file=23
  open(file, FILE="mxyz_solid.in", STATUS="old",action="read")

  do inn=1,nn
     read(file,*) xyz(inn,1:nsd), n_bc(inn), bc_type(1:n_bc(inn),inn)
  enddo
  close(file)
  return
end subroutine readx_solid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! added solid part number, mat_part, in the data for multiple types of materials in solid
!**************************************************************************************
subroutine readien_solid(solid_con,ne,nen,mat_part)
  implicit none

  integer,intent(in) :: ne,nen
  integer :: solid_con(ne,nen),mat_part(ne)

  integer :: file,ine


  file=21
  open(file, FILE="mien_solid.in", STATUS="old",action="read")

  do ine=1,ne
     read(file,*) solid_con(ine,1:nen), mat_part(ine)
  enddo
  close(file)

  return
end subroutine readien_solid

end module meshgen_solid
