module meshgen_solid
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readx_solid(xyz,nn,nsd)
  implicit none

  integer,intent(in) :: nn,nsd
  real*8 :: xyz(nn,nsd)

  integer :: idummy,inn,file


  file=23
  open(file, FILE="mxyz_solid.in", STATUS="old",action="read")

  do inn=1,nn
     read(file,*) idummy,xyz(inn,1:nsd),idummy
  enddo	
  close(file)
  write(*,*) "meshgen",nn

  return
end subroutine readx_solid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readien_solid(solid_con,ne,nen)
  implicit none

  integer,intent(in) :: ne,nen
  integer solid_con(ne,nen)

  integer :: file,ine,idummy


  file=21
  open(unit=file, FILE="mien_solid.in", STATUS="old",action="read")
  do ine=1,ne
     read(file,*) idummy,solid_con(ine,1:nen)
  enddo
  close(file)

  return
end subroutine readien_solid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readrng_solid(rngface,ne,neface)
  implicit none

  integer,intent(in) :: ne,neface
  integer rngface(ne,neface)
	!integer rngfaceieee(neface,ne/2+1)
!	integer lock,ierr,io,status(MPI_STATUS_SIZE)
	!character*4 ifp
  integer file,i,ieface,iec
	!integer offset,endset

  file=26
  open(file, FILE="mrng_solid.in", STATUS="old",action="read")
  do i=1,ne
     read(file,*) rngface(i,:)
  enddo

  do ieface=1,neface
     do iec=1,ne
        if(rngface(iec,ieface).lt.0) rngface(iec,ieface) = 0
     enddo
  enddo
	
  !mynrng = 0
  !do ieface=1,neface
  !   do iec=1,ne
  !      mynrng = max(mynrng, rngface(ieface,iec))
  !   end do
  !end do
  !nrng=mynrng
  close(file)

  return
end subroutine readrng_solid


end module meshgen_solid