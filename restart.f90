!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! modified form io11.f file to generate 
! ensight fluid field file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program restart !(restart_value,d)
!  use solid_variables, only:nn_solid
!  use fluid_variables, only: nn,ndf

  implicit none
!	real*8 :: d(ndf,nn),solid_coor_curr(3,nn_solid)

  real*8 :: d(4,4861),solid_coor_curr(3,19940),solid_vel(3,19940)
  integer :: klok

	integer:: nn_solid,restart_value,nn
	character*12::file_name
	character*80::dummy
	integer::i_file_unit,i,j

  integer :: in,n_sum
  character*5  fileroot
  character*12 name_file1
  character*12 name_file2
  character*15 name_file3
  character*15 name_file4

  integer,parameter :: ifileunit = 15

! read solid geometry file
! from
	nn_solid=19940
	nn=       4861
	restart_value=0
	i_file_unit=20

 200  format(a5   )
 201  format(a4,i1)
 202  format(a3,i2)
 203  format(a2,i3)
 204  format(a1,i4)
 205  format(   i5)
 klok = restart_value

  write(*,*) 'klok=',klok

  if (klok .eq. 0) then
     write(fileroot, 200) '00000'
  elseif (klok .lt. 10) then
     write(fileroot, 201) '0000',klok
  elseif (klok .lt. 100) then
     write(fileroot, 202) '000' ,klok
  elseif (klok .lt. 1000) then
     write(fileroot, 203) '00'  ,klok
  elseif (klok .lt. 10000) then
     write(fileroot, 204) '0'   ,klok
  elseif (klok .lt. 100000) then    
     write(fileroot, 205)   ''  ,klok
  else
     write(0,*) 'klok .ge. 100000: modify subroutine createfileroot'
     call exit(1)
  endif

  write(file_name,'(A7, A5)')  'fem.geo', fileroot

  open(i_file_unit, file=file_name, status='old',form='formatted')

  write(*,*) 'reading... ',file_name

  read(i_file_unit, *) dummy
  read(i_file_unit, *) dummy

  read(i_file_unit, *) dummy
  read(i_file_unit, *) dummy
  read(i_file_unit, *) dummy
  read(i_file_unit, '(I8)') n_sum

  read(i_file_unit,101) (j,solid_coor_curr(1,j),solid_coor_curr(2,j), &
       solid_coor_curr(3,j),j=1,nn_solid)

101	format(i8,3e12.5)
  write(*,*) solid_coor_curr(1,1),solid_coor_curr(2,1),solid_coor_curr(3,1)
  close(i_file_unit)

  write(name_file1,'(A7,  A5)')  'fem.vel', fileroot
!  write(name_file2,'(A7,  A5)')  'fem.pre', fileroot
!  write(name_file3,'(A10, A5)')  'fem.stress', fileroot
!  write(name_file4,'(A10, A5)')  'fem.strain', fileroot


 !...Write velocity output in ens_movie.vel*
  write(*,*) 'reading... ', name_file1

  open(ifileunit, file=name_file1, status='old', form='formatted')
  read(ifileunit, *) dummy
  read(ifileunit,110) (solid_vel(1,in),solid_vel(2,in),solid_vel(3,in),in=1,nn_solid), &
                       (d(1,in),d(2,in),d(3,in),in=1,nn)
  close(ifileunit)
  write(*,*) solid_vel(1,1),solid_vel(2,1),solid_vel(3,1),d(1,1),d(2,1),d(3,1)


!...Write pressure output in ens_movie.pre*
!  write(*,*) 'reading... ', name_file2
!  open(ifileunit, file=name_file2, form='formatted',status='old')
!  read(ifileunit, '(A)') 'structure and fluid field: pressure'  
!  read(ifileunit,110) (pave(in),in=1,nn_solid),(d(4,in),in=1,nn)
!  close(ifileunit)


110	format(6e12.5)

end program restart