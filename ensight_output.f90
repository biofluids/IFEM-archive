module ensight_output
  implicit none

contains

subroutine zfem_ensCase(dt, currentStep,ntsbout)
  implicit none

  real(8) :: dt
  integer :: currentStep,ntsbout

  integer :: numbers_of_step
  real(8) :: time_value(10000)

  character(len=12) :: file_name
  !character*13 mgeo_name
  character(len=12) :: pre_name
  character(len=12) :: vel_name
  character(len=12) :: FSI_name
  character(len=15) :: stress
  character(len=15) :: strain     

  integer :: ts 
  integer :: file_start_no, file_incre

  integer :: i,k

  write(*,*) 'generate ensight case file'

  ts = 1
  numbers_of_step = currentStep/ntsbout
  file_start_no = 0
  file_incre = ntsbout


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     write case file 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  open (unit=20, file = 'fem.case', status='unknown')
  write(20, '(A6)') 'FORMAT'
  write(20, '(A20)') 'type:   ensight'
  write(20, *) 
  write(20, '(A8)') 'GEOMETRY'

  file_name = 'fem.geo*****'
  !mgeo_name = 'fem.geo*****' 
  pre_name  = 'fem.pre*****'
  vel_name  = 'fem.vel*****'
  fsi_name  = 'fem.fsi*****'
  stress    = 'fem.stress*****'
  strain    = 'fem.strain*****'
 
 !5001 format(A9, 11x, i5, 5x, A13, 1x, A18)
 5002 format(A16, 1x, I2, 1x, A8, 1x, A)

 5100 format(A9, 21x, i5)
 5101 format(A16, 14x, i5)
 5102 format(A22, 8x, i5)
 5103 format(A19, 11x, i5)

      write(20, 5000) 'model:', ts, file_name, 'change_coords_only'
 5000 format(A6, 14x, i5, 5x, A12, 1x, A18) 

  write(20, *) 

  write(20, '(A8)') 'VARIABLE'
  write(20, 5002) 'scalar per node:', ts, 'pressure', pre_name 
  write(20, 5002) 'vector per node:', ts, 'velocity', vel_name
  write(20, 5002) 'vector per node:', ts, 'forceFSI', FSI_name

  write(20, 5003) ts, stress
  write(20, 5004) ts, strain 
 5003 format('tensor symm per node: ',I2, 1x,'stress ', A)
 5004 format('tensor symm per node: ',I2, 1x,'strain ', A)

  write(20, *)
  write(20, *)

  write(20, '(A4)') 'TIME'
  write(20, 5100) 'time set:', ts
  write(20, 5101) 'number of steps:', numbers_of_step+1
  write(20, 5102) 'filename start number:', file_start_no
  write(20, 5103) 'filename increment:', file_incre

  write(20, *)

  write(20, '(A12)') 'time values:'
  do i=1,numbers_of_step
     time_value(i)=dt*ntsbout*i
  enddo
  write(20,5110) 0.0,(time_value(k),k=1,numbers_of_step)
 5110 format(5f14.6 )

  close(20)

  return
end subroutine zfem_ensCase



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine generates Ensight6 geometry file
! Lucy Zhang
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine zfem_ensGeo(klok,ien,xn,solid_fem_con,solid_coor_curr)
  use solid_variables
  use fluid_variables
  implicit none

  integer :: klok
  integer :: ien(nen,ne)
  real(8)  :: xn(nsd,nn)
  integer,dimension(1:ne_solid,1:nen_solid) :: solid_fem_con   !...connectivity for solid FEM mesh
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_curr  !...node position current

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
     write(0,*) 'klok >= 100000: modify subroutine createfileroot'
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
101 format(i8,3e12.5)

 !...write fluids coordinates
  do i=1,nn
     if (nsd==3) write(i_file_unit,101) i+nn_solid,xn(1,i),xn(2,i),xn(3,i)
     if (nsd==2) write(i_file_unit,101) i+nn_solid,xn(1,i),xn(2,i),0.0
!     write(i_file_unit,101) i+nn_solid,xn(1:nsd,i)
  enddo


 

    !write fluids part element connectivity
  write(i_file_unit, *) 'part 1'
  write(i_file_unit, *) ' Fluid Model'
  if (nsd==3) then
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
  elseif (nsd==2) then
	  select case (nen)
	  case (3)
		 write(i_file_unit, *) ' tria3'  ! element type
		 write(i_file_unit, '(i8)')  ne   ! number of elements
		 do i=1, ne
			write(i_file_unit,'(4i8)') i, (ien(j,i)+nn_solid,j=1,nen) !element connectivity
		 enddo
	  case (4)
		 write(i_file_unit, *) 'quad4'    ! element type
		 write(i_file_unit, '(I8)')  ne   ! number of elements
		 do i=1, ne
			write(i_file_unit,'(9i8)') i, (ien(j,i)+nn_solid,j=1,nen) !element connectivity
		 enddo
	  end select
  endif


  close(i_file_unit)

  return
end subroutine zfem_ensGeo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! modified form io11.f file to generate 
! ensight fluid field file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine zfem_ensFluid(d,f_fluids,solid_force_FSI,solid_vel,solid_pave,solid_stress,solid_strain,klok)
  use solid_variables
  use fluid_variables, only: nn,ndf,nsd
  use run_variables, only: its
  implicit none

  real(8) :: d(ndf,nn)
  real(8) :: f_fluids(nsd,nn)

  !integer,dimension(1:ne_solid,1:nen_solid) :: solid_fem_con   !...connectivity for solid FEM mesh
  !integer,dimension(1:ne_solid,1:nsurface)  :: solid_surface   !...surface element faces

  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_force_FSI   !...fluid structure interaction force
  !real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_init   !...node position initial
  !real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_curr   !...node position current
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_vel         !...velocity
  !real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_prevel      !...velocity - previous timestep
  !real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_accel       !...acceleration


  real(8),dimension(nn_solid)   :: solid_pave  !...averaged solid pressure (from mixed formulation -> ???)

  real(8),dimension(6,nn_solid) :: solid_stress  !...solid stress (Voigt notation)
  real(8),dimension(6,nn_solid) :: solid_strain  !...solid strain (Voigt notation)
  integer :: klok
  
  real(8) :: fluid_stress(6,nn),fluid_strain(6,nn)              !...fluid stress and strain (not used) empty matrix needed for Ensight input format
  !real(8) :: solid_stress(6,nn_solid),solid_strain(6,nn_solid)    !...solid stress and strain

  !real(8) :: pave(nn_solid)

  integer :: in
  !integer :: ntem
  character(len=5 ) :: fileroot
  character(len=12) :: name_file1
  character(len=12) :: name_file2
  character(len=15) :: name_file3
  character(len=15) :: name_file4
  character(len=12) :: name_file5

  integer,parameter :: ifileunit = 15


  fluid_stress(1:6,1:nn)=0.0
  fluid_strain(1:6,1:nn)=0.0


  if (klok .eq. 0) then
     write(fileroot, '(a5)') '00000'
  elseif (klok .lt. 10) then
     write(fileroot, '(a4,i1)') '0000',klok
  elseif (klok .lt. 100) then
     write(fileroot, '(a3,i2)') '000' ,klok
  elseif (klok .lt. 1000) then
     write(fileroot, '(a2,i3)') '00'  ,klok
  elseif (klok .lt. 10000) then
     write(fileroot, '(a1,i4)') '0'   ,klok
  elseif (klok .lt. 100000) then    
     write(fileroot, '(i5)')   ''  ,klok
  else
     write(0,*) 'klok >= 100000: modify subroutine createfileroot'
     call exit(1)
  endif

  write(name_file1,'(A7,  A5)')  'fem.vel', fileroot
  write(name_file2,'(A7,  A5)')  'fem.pre', fileroot
  write(name_file3,'(A10, A5)')  'fem.stress', fileroot
  write(name_file4,'(A10, A5)')  'fem.strain', fileroot
  write(name_file5,'(A7,  A5)')  'fem.fsi', fileroot

!===========================================================================
! Output velocity, interaction force, pressure, stress, strain into Ensight format
!===========================================================================

!=========  3D =====================================
if (nsd==3) then
 !...Write velocity output in ens_movie.vel*
  write(*,*) 'writing... ', name_file1
  open(ifileunit, file=name_file1, form='formatted')
  write(ifileunit, '(A)')   'structure and fluid field: velocity vector'
  write(ifileunit,110) (solid_vel(1,in),solid_vel(2,in),solid_vel(3,in),in=1,nn_solid), &
                       (d(1,in),d(2,in),d(3,in),in=1,nn)
  close(ifileunit)

 !...Write Interaction force output in ens_movie.fsi*
  write(*,*) 'writing... ', name_file5
  open(ifileunit, file=name_file5, form='formatted')
  write(ifileunit, '(A)')   'structure and fluid field: force_FSI vector'
  write(ifileunit,110) (solid_force_FSI(1,in),solid_force_FSI(2,in),solid_force_FSI(3,in),in=1,nn_solid), &
                       (f_fluids(1,in),f_fluids(2,in),f_fluids(3,in),in=1,nn)
  close(ifileunit)

 !...Write pressure output in ens_movie.pre*
  write(*,*) 'writing... ', name_file2
  open(ifileunit, file=name_file2, form='formatted')
  write(ifileunit, '(A)') 'structure and fluid field: pressure'  
  write(ifileunit,110) (solid_pave(in),in=1,nn_solid),(d(4,in),in=1,nn)
  close(ifileunit)

 !...Write stress output in ens_movie.stress*
  write(*,*) 'writing... ', name_file3
  open(ifileunit, file=name_file3, form='formatted')
  write(ifileunit, '(A)') 'structure field: stress'  
  write(ifileunit,110) (solid_stress(1,in),solid_stress(2,in),solid_stress(3,in),                 &
                        solid_stress(4,in),solid_stress(5,in),solid_stress(6,in),in=1,nn_solid),  &
                       (fluid_stress(1,in),fluid_stress(2,in),fluid_stress(3,in),                 &
                        fluid_stress(4,in),fluid_stress(5,in),fluid_stress(6,in),in=1,nn)
  close(ifileunit)

 !...Write strain output in ens_movie.strain*
  write(*,*) 'writing... ', name_file4
  open(ifileunit, file=name_file4, form='formatted')
  write(ifileunit, '(A)') 'structure field: strain'  
  write(ifileunit,110) (solid_strain(1,in),solid_strain(2,in),solid_strain(3,in),                 &
                        solid_strain(4,in),solid_strain(5,in),solid_strain(6,in),in=1,nn_solid),  &
                       (fluid_strain(1,in),fluid_strain(2,in),fluid_strain(3,in),                 &
                        fluid_strain(4,in),fluid_strain(5,in),fluid_strain(6,in),in=1,nn)
  close(ifileunit)

110 format(6e12.5)

!===================  2D  =======================================
elseif (nsd==2) then
  !...Write velocity output in ens_movie.vel*
  write(*,*) 'writing... ', name_file1
  open(ifileunit, file=name_file1, form='formatted')
  write(ifileunit, '(A)')   'structure and fluid field: velocity vector'
  write(ifileunit,110) (solid_vel(1,in),solid_vel(2,in),0.0,in=1,nn_solid), &
                       (d(1,in),d(2,in),0.0,in=1,nn)
  close(ifileunit)

 !...Write Interaction force output in ens_movie.fsi*
  write(*,*) 'writing... ', name_file5
  open(ifileunit, file=name_file5, form='formatted')
  write(ifileunit, '(A)')   'structure and fluid field: force_FSI vector'
  write(ifileunit,110) (solid_force_FSI(1,in),solid_force_FSI(2,in),solid_force_FSI(3,in),in=1,nn_solid), &
                       (f_fluids(1,in),f_fluids(2,in),0.0,in=1,nn)
  close(ifileunit)

 !...Write pressure output in ens_movie.pre*
  write(*,*) 'writing... ', name_file2
  open(ifileunit, file=name_file2, form='formatted')
  write(ifileunit, '(A)') 'structure and fluid field: pressure'  
  write(ifileunit,110) (solid_pave(in),in=1,nn_solid),(d(3,in),in=1,nn)
  close(ifileunit)

 !...Write stress output in ens_movie.stress*
  write(*,*) 'writing... ', name_file3
  open(ifileunit, file=name_file3, form='formatted')
  write(ifileunit, '(A)') 'structure field: stress'  
  write(ifileunit,110) (solid_stress(1,in),solid_stress(2,in),solid_stress(3,in),                 &
                        0.0, 0.0, solid_stress(6,in),in=1,nn_solid),  &
                       (fluid_stress(1,in),fluid_stress(2,in),fluid_stress(3,in),                 &
                        0.0, 0.0,fluid_stress(6,in),in=1,nn)
  close(ifileunit)

 !...Write strain output in ens_movie.strain*
  write(*,*) 'writing... ', name_file4
  open(ifileunit, file=name_file4, form='formatted')
  write(ifileunit, '(A)') 'structure field: strain'  
  write(ifileunit,110) (solid_strain(1,in),solid_strain(2,in),solid_strain(3,in),                 &
                        0.0, 0.0,solid_strain(6,in),in=1,nn_solid),  &
                       (fluid_strain(1,in),fluid_strain(2,in),fluid_strain(3,in),                 &
                        0.0, 0.0,fluid_strain(6,in),in=1,nn)
  close(ifileunit)


endif

  return
end subroutine zfem_ensFluid


end module ensight_output
