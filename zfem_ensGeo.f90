!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine generates Ensight6 geometry file
! Lucy Zhang
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine zfem_ensGeo(klok,ien,xn)
  use solid_variables
  use r_common
  use fluid_variables
  implicit none

  integer :: klok
  integer :: ien(nen,ne)
  real*8  :: xn(nsd,nn)

  character(len =  5) :: fileroot
  character(len = 12) :: file_name

  integer,parameter :: i_file_unit = 15

!%%%%%%%%%%%%%%%%%%%%%%
! used for ensight
!%%%%%%%%%%%%%%%%%%%%%
  !integer sizeX, sizeY, sizeZ
  integer :: i, j



 200  format(a5   )
 201  format(a4,i1)
 202  format(a3,i2)
 203  format(a2,i3)
 204  format(a1,i4)
 205  format(   i5)

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
  write(6,*) 'writing... ', file_name 

  open(i_file_unit, file=file_name, form='formatted')
  write(i_file_unit, *) 'This is the ensight format geometry file'
  write(i_file_unit, *) 'This is the ensitht format geometry file'
      

  write(i_file_unit, *) 'node id given'
  write(i_file_unit, *) 'element id given'
  write(i_file_unit, *) 'coordinates'
  write(i_file_unit, '(I8)') nn_solid + nn

!--> node id, x, y, x
!-->
 !write structure coordinates
  do i=1, nn_solid
     write(i_file_unit, 101) i,solid_coor_curr(1,i),  &
	                           solid_coor_curr(2,i),  &
							   solid_coor_curr(3,i)
  enddo
101	format(i8,3e12.5)

 !...write fluids coordinates
  do i=1,nn
     write(i_file_unit,101) i+nn_solid,xn(1,i),xn(2,i),xn(3,i)
  enddo


 !...write structure part - connectivity
  write(i_file_unit, *) 'part 1'
  write(i_file_unit, *) ' Structure Model'
  if (nsd_solid == 0) then
     write(i_file_unit,'(a7)') '  point'
	 write(i_file_unit,'(i8)') nn_solid
	 do i=1,nn_solid
		write(i_file_unit,'(2i8)') i,i
	 enddo
  elseif (nsd_solid == 3) then
     select case (nis)
     case (8) 
        write(i_file_unit,'(a7)') '  hexa8'    ! element type
	    write(i_file_unit, '(i8)')  ne_solid   ! number of elements
        do i=1, ne_solid
	       write(i_file_unit,'(9i8)') i, (nea(i,j),j=1,nis) !element connectivity
	    enddo	  
     case (4)
	    write(i_file_unit,'(a7)') ' tetra4'    ! element type
	    write(i_file_unit, '(i8)')  ne_solid   ! number of elements
        do i=1, ne_solid
           write(i_file_unit,'(5i8)') i, (nea(i,j),j=1,nis) !element connectivity
		enddo
	 case default
	    write(*,*) "zfem_ens: no ensight output defined for nis = ",nis
		stop
	 end select
  endif

	!write fluids part element connectivity
  write(i_file_unit, *) 'part 2'
  write(i_file_unit, *) ' Fluid Model'
  select case (nen)
  case (4)
     write(i_file_unit, *) ' tetra4'  ! element type
     write(i_file_unit, '(i8)')  ne   ! number of elements
     do i=1, ne
        write(i_file_unit,'(5i8)') i, (ien(j,i)+nn_solid,j=1,nen) !element connectivity
	 enddo
  case (8)
     write(i_file_unit, *) 'hexa8'    ! element type
     write(i_file_unit, '(I8)')  ne   ! number of elements
     do i=1, ne
        write(i_file_unit,'(9i8)') i, (ien(j,i)+nn_solid,j=1,nen) !element connectivity
     enddo
  end select


  close(i_file_unit)

  return
end subroutine zfem_ensGeo
