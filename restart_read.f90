!===========================================
! Read data from ensight output for restart
! Lucy Zhang, 5/27/03
!===========================================

!================================================
! Read solid geometry file from fem.geo*****
!================================================
  subroutine read_x (restart_value,solid_coor_curr)

  use solid_variables, only:nn_solid,nsd_solid
!  use fluid_variables, only: nn,ndf

  implicit none
  integer :: klok

  character*12::file_name
  character*80::dummy
  integer::restart_value, j
  real* 8 solid_coor_curr(nsd_solid,nn_solid)

  integer :: n_sum
  character*5  fileroot

  integer,parameter :: ifileunit = 15

 200  format(a5   )
 201  format(a4,i1)
 202  format(a3,i2)
 203  format(a2,i3)
 204  format(a1,i4)
 205  format(   i5)
 klok = restart_value

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

  open(ifileunit, file=file_name, status='old',form='formatted')

  write(*,*) 'reading... ',file_name

  read(ifileunit, *) dummy
  read(ifileunit, *) dummy

  read(ifileunit, *) dummy
  read(ifileunit, *) dummy
  read(ifileunit, *) dummy
  read(ifileunit, '(I8)') n_sum

  read(ifileunit,101) (j,solid_coor_curr(1,j),solid_coor_curr(2,j), &
       solid_coor_curr(3,j),j=1,nn_solid)

101	format(i8,3e12.5)
!  write(*,*) solid_coor_curr(1,1),solid_coor_curr(2,1),solid_coor_curr(3,1)
  close(ifileunit)

  return
  end subroutine read_x


!=====================================================
! read solid and fluid velocity file fem.vel*****
!=====================================================
  subroutine read_vel (restart_value,solid_vel,d,solid_accel)

  use solid_variables, only:nn_solid
  use fluid_variables, only: nn,ndf,d

  implicit none
  real*8 :: d(ndf,nn),solid_vel(3,nn_solid),solid_accel(3,nn_solid)
  integer :: klok
  character*80::dummy
  integer::restart_value

  integer :: in,i
  character*5  fileroot
  character*12 name_file1

  integer,parameter :: ifileunit = 15

 200  format(a5   )
 201  format(a4,i1)
 202  format(a3,i2)
 203  format(a2,i3)
 204  format(a1,i4)
 205  format(   i5)
 klok = restart_value

  write(*,*) 'READING RESULTS FROM TIME STEP=',klok

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

  write(name_file1,'(A7,  A5)')  'fem.vel', fileroot

 !...Write velocity output in ens_movie.vel*
  write(*,*) 'reading... ', name_file1

  open(ifileunit, file=name_file1, status='old', form='formatted')
  read(ifileunit, *) dummy
  read(ifileunit,110) (solid_vel(1,in),solid_vel(2,in),solid_vel(3,in),in=1,nn_solid), &
                       (d(1,in),d(2,in),d(3,in),in=1,nn)
  close(ifileunit)

  open(222, file='accel.dat', status='old', form='formatted')
  read(222,110) (solid_accel(1,i),solid_accel(2,i),solid_accel(3,i),i=1,nn_solid)
  close(222)

110	format(6e12.5)
 return
 end subroutine read_vel