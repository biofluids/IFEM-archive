!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module: output_gmsh.f90
!
!  Axel Gerstenberger, NWU, Oct. 16 2003
!
!  provides subroutines to write output in the Gmsh post-processing format
!  for more information on the output format see the Gmsh documentation
!
module output_gmsh
  implicit none
  save


contains

subroutine write_1d_scalar(nn,nsd,x,scal,str_length,str,file_str_length,file_str)
  implicit none
  
  integer,intent(in) :: nn,nsd
  real(8),intent(in) :: x(1:nsd,1:nn)
  real(8),intent(in) :: scal(1:nn)
  integer,intent(in)             :: str_length
  character(len=str_length)      :: str
  integer,intent(in)             :: file_str_length
  character(len=file_str_length) :: file_str

  integer,parameter :: file = 21
  integer :: in

  open(file, file=file_str, form="formatted", status="replace")

  write(file,*) 'View "',str,'" {'
  do in = 1,nn
     write(unit=file,fmt='("SP(",F12.8,",",F12.8,",",F12.8,"){",F16.8,"};")') x(1,in),x(2,in),x(3,in),scal(in)
  enddo
  write(file,*) "};"

  close(file)

end subroutine write_1d_scalar


subroutine write_2d_scalar(nn,ne,nen,nsd,ien,x,scal,str_length,str,file_str_length,file_str)
  implicit none
  
  integer,intent(in) :: nn,ne,nen,nsd
  integer,intent(in) :: ien(1:nen,1:ne)
  real(8),intent(in) :: x(1:nsd,1:nn)
  real(8),intent(in) :: scal(1:nn)
  integer,intent(in)             :: str_length
  character(len=str_length)      :: str
  integer,intent(in)             :: file_str_length
  character(len=file_str_length) :: file_str

  integer,parameter :: file = 21

  integer :: ine
  

  open(file, file=file_str, form="formatted", status="replace")
  write(file,*) 'View "',str,'" {'
  do ine = 1,ne
     write(unit=file,fmt='("ST(",E16.8,",",E16.8,",",E16.8,"," ,   &
                                 E16.8,",",E16.8,",",E16.8,"," ,   &
                                 E16.8,",",E16.8,",",E16.8,"){",   &
                                 E16.8,",",E16.8,",",E16.8,"};")')    &
     x(1,ien(1,ine)),x(2,ien(1,ine)),x(3,ien(1,ine)),     &
     x(1,ien(2,ine)),x(2,ien(2,ine)),x(3,ien(2,ine)),     &
     x(1,ien(3,ine)),x(2,ien(3,ine)),x(3,ien(3,ine)),     &
     scal(ien(1,ine)),scal(ien(2,ine)),scal(ien(3,ine))
  enddo
  
  write(file,*) "};"

  close(file)

end subroutine write_2d_scalar

subroutine write_2d_vector(nn,ne,nen,nsd,ien,x,vect,str_length,str,file_str_length,file_str)
  implicit none
  
  integer,intent(in) :: nn,ne,nen,nsd
  integer,intent(in) :: ien(1:nen,1:ne)
  real(8),intent(in) :: x(1:nsd,1:nn)
  real(8),intent(in) :: vect(1:nsd,1:nn)
  integer,intent(in)             :: str_length
  character(len=str_length)      :: str
  integer,intent(in)             :: file_str_length
  character(len=file_str_length) :: file_str

  integer,parameter :: file = 21

  integer :: ine
  

  open(file, file=file_str, form="formatted", status="replace")
  write(file,*) 'View "',str,'" {'
  do ine = 1,ne
     write(unit=file,fmt='("VT(",F12.8,",",F12.8,",",F12.8,"," ,   &
                                 F12.8,",",F12.8,",",F12.8,"," ,   &
                                 F12.8,",",F12.8,",",F12.8,"){",   &
                                 F12.8,",",F12.8,",",F12.8,"," ,   &
                                 F12.8,",",F12.8,",",F12.8,"," ,   &
                                 F12.8,",",F12.8,",",F12.8,"};")')    &
     x(1,ien(1,ine)),x(2,ien(1,ine)),x(3,ien(1,ine)),     &
     x(1,ien(2,ine)),x(2,ien(2,ine)),x(3,ien(2,ine)),     &
     x(1,ien(3,ine)),x(2,ien(3,ine)),x(3,ien(3,ine)),     &
     vect(1,ien(1,ine)),vect(2,ien(1,ine)),vect(3,ien(1,ine)),  &
     vect(1,ien(2,ine)),vect(2,ien(2,ine)),vect(3,ien(2,ine)),  &
     vect(1,ien(3,ine)),vect(2,ien(3,ine)),vect(3,ien(3,ine))
  enddo
  
  write(file,*) "};"

  close(file)

end subroutine write_2d_vector



subroutine write_3d_scalar(nn,ne,nen,nsd,ien,x,scal,str_length,str,file_str_length,file_str)
  implicit none
  
  integer,intent(in) :: nn,ne,nen,nsd
  integer,intent(in) :: ien(1:nen,1:ne)
  real(8),intent(in) :: x(1:nsd,1:nn)
  real(8),intent(in) :: scal(1:nn)
  integer,intent(in)             :: str_length
  character(len=str_length)      :: str
  integer,intent(in)             :: file_str_length
  character(len=file_str_length) :: file_str

  integer,parameter :: file = 21

  integer :: ine
  

  open(file, file=file_str, form="formatted", status="replace")
  write(file,*) 'View "',str,'" {'
  do ine = 1,ne
     write(unit=file,fmt='("SS(",F12.8,",",F12.8,",",F12.8,"," ,   &
                                 F12.8,",",F12.8,",",F12.8,"," ,   &
                                 F12.8,",",F12.8,",",F12.8,"," ,   &
                                 F12.8,",",F12.8,",",F12.8,"){",   &
                                 F12.8,",",F12.8,",",F12.8,",",F12.8,"};")')    &
     x(1,ien(1,ine)),x(2,ien(1,ine)),x(3,ien(1,ine)),     &
     x(1,ien(2,ine)),x(2,ien(2,ine)),x(3,ien(2,ine)),     &
     x(1,ien(3,ine)),x(2,ien(3,ine)),x(3,ien(3,ine)),     &
     x(1,ien(4,ine)),x(2,ien(4,ine)),x(3,ien(4,ine)),     &
     scal(ien(1,ine)),scal(ien(2,ine)),scal(ien(3,ine)),scal(ien(4,ine))
  enddo
  
  write(file,*) "};"

  close(file)

end subroutine write_3d_scalar

end module output_gmsh