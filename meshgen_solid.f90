module meshgen_solid
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readx_solid(xyz,nn,nsd)
  use mpi_variables
  implicit none

  integer,intent(in) :: nn,nsd
  real(8) :: xyz(nn,nsd)
  integer :: inn,file

  file=23
  open(file, FILE="mxyz_solid.in", STATUS="old",action="read")
if (nsd == 2) then
  do inn=1,nn
     read(file,100) xyz(inn,1:nsd)
  enddo
end if

if (nsd == 3) then
  do inn=1,nn
     read(file,101) xyz(inn,1:nsd)
if(myid==0)write(*,*)'nn=',inn,xyz(inn,1:nsd)
  enddo
end if



100 format (D14.10, D14.10)
101 format (D14.10, D14.10, D14.10)

  close(file)

  return
end subroutine readx_solid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readien_solid(solid_con,ne,nen,mtype)
  implicit none

  integer,intent(in) :: ne,nen
  integer :: solid_con(ne,nen)
  integer :: mtype(ne)

  integer :: file,ine

  file=21
  open(file, FILE="mien_solid.in", STATUS="old",action="read")
if (nen == 3) then
  do ine=1,ne
     read(file,100) solid_con(ine,1:nen), mtype(ine)
  enddo
end if

if (nen == 4) then
  do ine=1,ne
     read(file,200) solid_con(ine,1:nen), mtype(ine)
  enddo
end if

if (nen == 8) then
  do ine=1,ne
     read(file,300) solid_con(ine,1:nen), mtype(ine)
  enddo
end if


100 format (I8,I8,I8,I8)
200 format (I8,I8,I8,I8,I8)
300 format(I8,I8,I8,I8,I8,I8,I8,I8,I8)
  close(file)

  return
end subroutine readien_solid
!---------------------------------------------------------------
subroutine read_sfcon(sfcon,node_sfcon)
  implicit none

  integer,intent(in) :: node_sfcon
  integer :: sfcon(node_sfcon)
  integer :: file,i

  file=23
  open(file, FILE="sfcon.in", STATUS="old",action="read")

  do i=1,node_sfcon
     read(file,100) sfcon(i)
  enddo
100 format (I8)
  close(file)

  return
end subroutine read_sfcon

end module meshgen_solid
