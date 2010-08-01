module meshgen_fluid
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readx(xyz)
  use fluid_variables, only: nsd,nn
  implicit none

  real(8) :: xyz(nsd,nn)
  integer :: file,i

  file=21
  open(file, FILE="mxyz.in", STATUS="old")

  do i=1,nn
     read(file,100) xyz(1:nsd,i)
  enddo
  close(file)
100 format(D14.10,D14.10)
  return
end subroutine readx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readien(ien)
  use fluid_variables, only: nen,ne
  implicit none

  integer :: ien(nen,ne)
  integer :: file,i

  file=21
  open(file, FILE="mien.in", STATUS="old")
if (nen==3) then
  do i=1,ne
     read(file,100) ien(:,i)
  enddo
end if

if (nen==4) then
  do i=1,ne
     read(file,101) ien(:,i)
  enddo
end if
  close(file)
100 format(I8,I8,I8)
101 format(I8,I8,I8,I8)
  return
end subroutine readien

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readrng(rngface)
  use fluid_variables, only: ne,neface,nrng,nen
  implicit none

  integer :: rngface(neface,ne)
  integer :: file,i,ieface,iec,mynrng

  file=21
  open(file, FILE="mrng.in", STATUS="old")
if (nen==3) then
  do i=1,ne
     read(file,100) rngface(:,i)
  enddo
end if

if (nen==4) then
  do i=1,ne
     read(file,101) rngface(:,i)
  enddo
end if
  do ieface=1,neface
     do iec=1,ne
        if(rngface(ieface,iec).lt.0) rngface(ieface,iec) = 0
     enddo
  enddo

  mynrng = 0
  do ieface=1,neface
     do iec=1,ne
        mynrng = max(mynrng, rngface(ieface,iec))
     end do
  end do
  nrng=mynrng
  close(file)
100 format(I8,I8,I8)
101 format(I8,I8,I8,I8)
  return
end subroutine readrng


end module meshgen_fluid
