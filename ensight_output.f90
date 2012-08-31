module ensight_output
  implicit none

contains

subroutine zfem_ensCase(dt, currentStep,ntsbout)
  implicit none

  real(8) :: dt
  integer :: currentStep,ntsbout

  integer :: numbers_of_step
  real(8) :: time_value(100000)

  character(len=13) :: file_name
  character(len=13) :: pre_name
  character(len=13) :: vel_name
  character(len=13) :: FSI_name
  character(len=16) :: stress
  character(len=16) :: strain     
  character(len=13) :: I_fluid
  character(len=13) :: I_solid
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

  file_name = 'fem.geo******'
  pre_name  = 'fem.pre******'
  vel_name  = 'fem.vel******'
  fsi_name  = 'fem.fsi******'
  stress    = 'fem.stress******'
  strain    = 'fem.strain******'
!  Indicator = 'fem.ind******'
  I_fluid = 'fem.inf******'
  I_solid = 'fem.ins******'


 5002 format(A16, 1x, I2, 1x, A8, 1x, A)
 5100 format(A9, 21x, i5)
 5101 format(A16, 14x, i5)
 5102 format(A22, 8x, i5)
 5103 format(A19, 11x, i5)
 5999 format(A16, 1x, I2, 1x, A7, 1x, A)
      write(20, 5000) 'model:', ts, file_name!, 'change_coords_only'
 5000 format(A6, 14x, i5, 5x, A13)!, 1x, A18) 

  write(20, *) 

  write(20, '(A8)') 'VARIABLE'
  write(20, 5002) 'scalar per node:', ts, 'pressure', pre_name 
  write(20, 5002) 'vector per node:', ts, 'velocity', vel_name
  write(20, 5002) 'vector per node:', ts, 'forceFSI', FSI_name
  write(20, 5999) 'scalar per node:', ts, 'I_solid', I_solid
  write(20, 5999) 'scalar per node:', ts, 'I_fluid', I_fluid

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
 5110 format(5f14.10 )

  close(20)

  return
end subroutine zfem_ensCase



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine generates Ensight6 geometry file
! Lucy Zhang
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine zfem_ensGeo(klok,ien,xn,solid_fem_con,solid_coor_curr,mtype,necover,nebase,x_inter)
  use solid_variables
  use fluid_variables
  use interface_variables, only:nn_inter,maxmatrix
  implicit none

  integer :: klok
  integer :: ien(nen,ne)
  real(8)  :: xn(nsd,nn)
  integer,dimension(1:ne_solid,1:nen_solid) :: solid_fem_con   !...connectivity for solid FEM mesh
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_curr  !...node position current
  integer :: mtype(ne)       ! Materail type index
  integer :: necover
  integer :: nebase
  real(8) :: x_inter(nsd,maxmatrix)
  character(len =  7) :: fileroot
  character(len = 14) :: file_name

  integer,parameter :: i_file_unit = 20

!%%%%%%%%%%%%%%%%%%%%%%
! used for ensight
!%%%%%%%%%%%%%%%%%%%%%
  integer :: i, j

 200  format(a6   )
 201  format(a5,i1)
 202  format(a4,i2)
 203  format(a3,i3)
 204  format(a2,i4)
 205  format(a1,i5)
 206  format(   i6)

  if (klok .eq. 0) then
     write(fileroot, 200) '000000'
  elseif (klok .lt. 10) then
     write(fileroot, 201) '00000',klok
  elseif (klok .lt. 100) then
     write(fileroot, 202) '0000' ,klok
  elseif (klok .lt. 1000) then
     write(fileroot, 203) '000'  ,klok
  elseif (klok .lt. 10000) then
     write(fileroot, 204) '00'   ,klok
  elseif (klok .lt. 100000) then
     write(fileroot, 205) '0'   ,klok
  elseif (klok .lt. 1000000) then    
     write(fileroot, 206)   ''  ,klok
  else
     write(0,*) 'klok >= 1000000: modify subroutine createfileroot'
     call exit(1)
  endif


  write(file_name,'(A7, A6)')  'fem.geo', fileroot
  write(6,*) 'writing... ', file_name 

  open(i_file_unit, file=file_name, form='formatted')
  write(i_file_unit, *) 'This is the ensight format geometry file'
  write(i_file_unit, *) 'This is the ensitht format geometry file'
      

  write(i_file_unit, *) 'node id given'
  write(i_file_unit, *) 'element id given'
  write(i_file_unit, *) 'coordinates'
  write(i_file_unit, '(I8)') nn_solid + nn + nn_inter
  
!--> node id, x, y, x
!-->
 !write structure coordinates
  do i=1, nn_solid
     if (nsd_solid==3) write(i_file_unit, 101) i,solid_coor_curr(1,i),  &
                               solid_coor_curr(2,i),  &
                               solid_coor_curr(3,i)
	 if (nsd_solid==2) write(i_file_unit, 101) i,solid_coor_curr(1,i),  &
                               solid_coor_curr(2,i), 0.0
  enddo
101 format(i8,3e12.5)

 !...write fluids coordinates
  do i=1,nn
     if (nsd==3) write(i_file_unit,101) i+nn_solid,xn(1,i),xn(2,i),xn(3,i)
     if (nsd==2) write(i_file_unit,101) i+nn_solid,xn(1,i),xn(2,i),0.0
  enddo

 !...write interface coordinates after regeneration
  do i=1,nn_inter
     if (nsd==3) write(i_file_unit,101) i+nn_solid+nn,x_inter(1,i),x_inter(2,i),x_inter(3,i)
     if (nsd==2) write(i_file_unit,101) i+nn_solid+nn,x_inter(1,i),x_inter(2,i),0.0
  end do



 !...write cover part - connectivity
  write(i_file_unit, *) 'part 1'
  write(i_file_unit, *) ' Cover '
  if (nsd_solid == 0) then
     write(i_file_unit,'(a7)') '  point'
     write(i_file_unit,'(i8)') nn_solid
     do i=1,nn_solid
        write(i_file_unit,'(2i8)') i,i
     enddo
  elseif (nsd_solid == 3) then
     select case (nen_solid)
     case (8) 
        write(i_file_unit,'(a7)') '  hexa8'    ! element type
        write(i_file_unit, '(i8)')  necover   ! number of elements
        do i=1, ne_solid
	   if (mtype(i) == 1) then
           write(i_file_unit,'(9i8)') i, (solid_fem_con(i,j),j=1,nen_solid) !element connectivity
	   end if
        enddo     
     case (4)
        write(i_file_unit,'(a7)') ' tetra4'    ! element type
        write(i_file_unit, '(i8)')  necover   ! number of elements
        do i=1, ne_solid
	   if (mtype(i) == 1) then
           write(i_file_unit,'(5i8)') i, (solid_fem_con(i,j),j=1,nen_solid) !element connectivity
	   end if
        enddo
     case default
        write(*,*) "zfem_ens: no ensight output defined for nen_solid = ",nen_solid
        stop
     end select
  elseif (nsd_solid == 2) then
	  select case (nen_solid)
	  case (3)
		 write(i_file_unit, *) ' tria3'  ! element type
		 write(i_file_unit, '(i8)')  necover   ! number of elements
		 do i=1, ne_solid
                    if (mtype(i)==1) then
			write(i_file_unit,'(4i8)') i, (solid_fem_con(i,j),j=1,nen_solid) !element connectivity
                    endif
		 enddo
	  case (4)
		 write(i_file_unit, *) 'quad4'    ! element type
		 write(i_file_unit, '(I8)')  necover   ! number of elements
		 do i=1, ne_solid
                   if (mtype(i)==1) then
			write(i_file_unit,'(5i8)') i, (solid_fem_con(i,j),j=1,nen_solid) !element connectivity
		   endif 
                enddo
	  end select
  endif
! Output solid base part
 write(i_file_unit, *) 'part 2'
  write(i_file_unit, *) ' Base '
if (nsd_solid == 0) then
     write(i_file_unit,'(a7)') '  point'
     write(i_file_unit,'(i8)') nn_solid
     do i=1,nn_solid
        write(i_file_unit,'(2i8)') i,i
     enddo
  elseif (nsd_solid == 3) then
     select case (nen_solid)
     case (8) 
        write(i_file_unit,'(a7)') '  hexa8'    ! element type
        write(i_file_unit, '(i8)')  nebase   ! number of elements
        do i=1, ne_solid
	   if (mtype(i) == 2) then
           write(i_file_unit,'(9i8)') i, (solid_fem_con(i,j),j=1,nen_solid) !element connectivity
	   end if
        enddo
     case (4)
        write(i_file_unit,'(a7)') ' tetra4'    ! element type
        write(i_file_unit, '(i8)')  nebase   ! number of elements
        do i=1, ne_solid
	   if (mtype(i) == 2) then
           write(i_file_unit,'(5i8)') i, (solid_fem_con(i,j),j=1,nen_solid) !element connectivity
	   end if
        enddo
     case default
        write(*,*) "zfem_ens: no ensight output defined for nen_solid = ",nen_solid
        stop
     end select
  elseif (nsd_solid == 2) then
          select case (nen_solid)
          case (3)
                 write(i_file_unit, *) ' tria3'  ! element type
                 write(i_file_unit, '(i8)')  nebase   ! number of elements
                 do i=1, ne_solid
                       if (mtype(i)==2) then
                        write(i_file_unit,'(4i8)') i, (solid_fem_con(i,j),j=1,nen_solid) !element connectivity
                       endif
                 enddo
          case (4)
                 write(i_file_unit, *) 'quad4'    ! element type
                 write(i_file_unit, '(I8)')  nebase   ! number of elements
                 do i=1, ne_solid
                    if (mtype(i)==2) then
                       write(i_file_unit,'(5i8)') i, (solid_fem_con(i,j),j=1,nen_solid) !element connectivity
                    endif 
                enddo
          end select
  endif

    !write fluids part element connectivity
  write(i_file_unit, *) 'part 3'
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
			write(i_file_unit,'(5i8)') i, (ien(j,i)+nn_solid,j=1,nen) !element connectivity
		 enddo
	  end select
  endif

  write(i_file_unit,*) 'part 4'
  write(i_file_unit,*) ' Interface Model'

  write(i_file_unit,'(a7)') '  point'
  write(i_file_unit,'(i8)') nn_inter
  do i=1,nn_inter
     write(i_file_unit,'(2i8)')i,i+nn_solid+nn
  end do



  close(i_file_unit)

  return
end subroutine zfem_ensGeo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! modified form io11.f file to generate 
! ensight fluid field file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine zfem_ensFluid(d,f_fluids,solid_force_FSI,solid_vel,solid_pave,solid_stress,solid_strain,klok,f_stress, &
         I_solid,I_fluid,vel_inter)

  use solid_variables
  use fluid_variables, only: nn,ndf,nsd,vis_liq
  use interface_variables, only:nn_inter,maxmatrix
  use run_variables, only: its
  implicit none

  real(8) :: d(ndf,nn)
  real(8) :: f_fluids(nsd,nn)
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_force_FSI   !...fluid structure interaction force
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_vel         !...velocity
  real(8),dimension (nsd,nsd,nn) :: f_stress
  real(8),dimension(nn_solid)   :: solid_pave  !...averaged solid pressure (from mixed formulation -> ???)
  real(8),dimension(1:nsd_solid*2,nn_solid) :: solid_stress  !...solid stress (Voigt notation)
  real(8),dimension(1:nsd_solid*2,nn_solid) :: solid_strain  !...solid strain (Voigt notation)
  real(8) I_solid(nn),I_fluid(nn),vel_inter(nsd,maxmatrix)
  integer :: klok
  real(8) :: fluid_stress(1:nsd*2,nn),fluid_strain(1:nsd*2,nn)              !...fluid stress and strain (not used) empty matrix needed for Ensight input format
  integer :: in,i
  character(len=7 ) :: fileroot
  character(len=14) :: name_file1
  character(len=14) :: name_file2
  character(len=17) :: name_file3
  character(len=17) :: name_file4
  character(len=14) :: name_file5
  character(len=14) :: name_file6
  character(len=14) :: name_file7
  integer,parameter :: ifileunit = 16

  ! Fluid stress and strain rate, Voigt notation
  fluid_stress(1:nsd*2,1:nn)=0.0
  fluid_strain(1:nsd*2,1:nn)=0.0

  if (nsd.eq.2) then
	! 2 Dim
	! Ensight Tensor Order: 1->11,2->22,3->33,4->12,5->13,6->23
	fluid_stress(1,1:nn)=f_stress(1,1,1:nn)		
	fluid_stress(2,1:nn)=f_stress(2,2,1:nn)		
	fluid_stress(4,1:nn)=f_stress(1,2,1:nn)		
	fluid_strain(1,1:nn)=f_stress(1,1,1:nn)/(2*vis_liq)			
	fluid_strain(2,1:nn)=f_stress(2,2,1:nn)/(2*vis_liq)
	fluid_strain(4,1:nn)=f_stress(1,2,1:nn)/(2*vis_liq)


  else
	! 3 Dim
	! Ensight Tensor Order: 1->11,2->22,3->33,4->12,5->13,6->23
	fluid_stress(1,1:nn)=f_stress(1,1,1:nn)
	fluid_stress(2,1:nn)=f_stress(2,2,1:nn)
	fluid_stress(3,1:nn)=f_stress(3,3,1:nn)
	fluid_stress(4,1:nn)=f_stress(1,2,1:nn)
	fluid_stress(5,1:nn)=f_stress(1,3,1:nn)
	fluid_stress(6,1:nn)=f_stress(2,3,1:nn)

	fluid_strain(1,1:nn)=f_stress(1,1,1:nn)/(2*vis_liq)
	fluid_strain(2,1:nn)=f_stress(2,2,1:nn)/(2*vis_liq)
	fluid_strain(3,1:nn)=f_stress(3,3,1:nn)/(2*vis_liq)
	fluid_strain(4,1:nn)=f_stress(1,2,1:nn)/(2*vis_liq)
	fluid_strain(5,1:nn)=f_stress(1,3,1:nn)/(2*vis_liq)
	fluid_strain(6,1:nn)=f_stress(2,3,1:nn)/(2*vis_liq)
  endif

  if (klok .eq. 0) then
     write(fileroot, '(a6)') '000000'
  elseif (klok .lt. 10) then
     write(fileroot, '(a5,i1)') '00000',klok
  elseif (klok .lt. 100) then
     write(fileroot, '(a4,i2)') '0000' ,klok
  elseif (klok .lt. 1000) then
     write(fileroot, '(a3,i3)') '000'  ,klok
  elseif (klok .lt. 10000) then
     write(fileroot, '(a2,i4)') '00'   ,klok
  elseif (klok .lt. 100000) then
     write(fileroot, '(a1,i5)') '0'   ,klok
  elseif (klok .lt. 1000000) then    
     write(fileroot, '(i6)')   ''  ,klok
  else
     write(0,*) 'klok >= 1000000: modify subroutine createfileroot'
     call exit(1)
  endif


  write(name_file1,'(A7,  A6)')  'fem.vel', fileroot
  write(name_file2,'(A7,  A6)')  'fem.pre', fileroot
  write(name_file3,'(A10, A6)')  'fem.stress', fileroot
  write(name_file4,'(A10, A6)')  'fem.strain', fileroot
  write(name_file5,'(A7,  A6)')  'fem.fsi', fileroot
  write(name_file6,'(A7,  A6)')  'fem.ins', fileroot
  write(name_file7,'(A7,  A6)')  'fem.inf', fileroot
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
                       (d(1,in),d(2,in),d(3,in),in=1,nn), &
                       (vel_inter(1,in),vel_inter(2,in),vel_inter(3,in),in=1,nn_inter)
  close(ifileunit)

 !...Write Interaction force output in ens_movie.fsi*
  write(*,*) 'writing... ', name_file5
  open(ifileunit, file=name_file5, form='formatted')
  write(ifileunit, '(A)')   'structure and fluid field: force_FSI vector'
  write(ifileunit,110) (solid_force_FSI(1,in),solid_force_FSI(2,in),solid_force_FSI(3,in),in=1,nn_solid), &
                       (f_fluids(1,in),f_fluids(2,in),f_fluids(3,in),in=1,nn), &
                       (0.0,0.0,0.0,in=1,nn_inter)
  close(ifileunit)

 !...Write pressure output in ens_movie.pre*
  write(*,*) 'writing... ', name_file2
  open(ifileunit, file=name_file2, form='formatted')
  write(ifileunit, '(A)') 'structure and fluid field: pressure'  
  write(ifileunit,110) (solid_pave(in),in=1,nn_solid),(d(4,in),in=1,nn),(0.0,in=1,nn_inter)
  close(ifileunit)

 !...Write stress output in ens_movie.stress*
  write(*,*) 'writing... ', name_file3
  open(ifileunit, file=name_file3, form='formatted')
  write(ifileunit, '(A)') 'structure field: stress'  
  ! Solid solver tensor order: 1->11,2->22,3->33,4->23,5->13,6->12
  ! Ensight tensor order: 1->11,2->22,3->33,4->12,5->13,6->23
  write(ifileunit,110) (solid_stress(1,in),solid_stress(2,in),solid_stress(3,in),                 &
                        solid_stress(6,in),solid_stress(5,in),solid_stress(4,in),in=1,nn_solid),  &
                       (fluid_stress(1,in),fluid_stress(2,in),fluid_stress(3,in),                 &
                        fluid_stress(4,in),fluid_stress(5,in),fluid_stress(6,in),in=1,nn), &
	               (0.0,0.0,0.0,0.0,0.0,0.0,in=1,nn_inter)
  close(ifileunit)

 !...Write strain output in ens_movie.strain*
  write(*,*) 'writing... ', name_file4
  open(ifileunit, file=name_file4, form='formatted')
  write(ifileunit, '(A)') 'structure field: strain'  
  ! Solid solver tensor order: 1->11,2->22,3->33,4->23,5->13,6->12
  ! Ensight tensor order: 1->11,2->22,3->33,4->12,5->13,6->23
  write(ifileunit,110) (solid_strain(1,in),solid_strain(2,in),solid_strain(3,in),                 &
                        solid_strain(6,in),solid_strain(5,in),solid_strain(4,in),in=1,nn_solid),  &
                       (fluid_strain(1,in),fluid_strain(2,in),fluid_strain(3,in),                 &
                        fluid_strain(4,in),fluid_strain(5,in),fluid_strain(6,in),in=1,nn), &
		       (0.0,0.0,0.0,0.0,0.0,0.0, in=1,nn_inter)
  close(ifileunit)

110 format(6e12.5)

  write(*,*) 'writing... ', name_file6
  open(ifileunit, file=name_file6, form='formatted')
  write(ifileunit, '(A)') 'structure and fluid field: I_solid'
  write(ifileunit,110) (0.0,in=1,nn_solid),(I_solid(in),in=1,nn),(0.5,in=1,nn_inter)
  close(ifileunit)

  write(*,*) 'writing... ', name_file7
  open(ifileunit, file=name_file7, form='formatted')
  write(ifileunit, '(A)') 'structure and fluid field: I_fluid'
  write(ifileunit,110) (0.0,in=1,nn_solid),(I_fluid(in),in=1,nn),(0.5,in=1,nn_inter)
  close(ifileunit)
!===================  2D  =======================================
elseif (nsd==2) then
  !...Write velocity output in ens_movie.vel*
  write(*,*) 'writing... ', name_file1
  open(ifileunit, file=name_file1, form='formatted')
  write(ifileunit, '(A)')   'structure and fluid field: velocity vector'
  write(ifileunit,110) (solid_vel(1,in),solid_vel(2,in),0.0,in=1,nn_solid), &
                       (d(1,in),d(2,in),0.0,in=1,nn), &
                       (vel_inter(1,in),vel_inter(2,in),0.0,in=1,nn_inter)
  close(ifileunit)

 !...Write Interaction force output in ens_movie.fsi*
  write(*,*) 'writing... ', name_file5
  open(ifileunit, file=name_file5, form='formatted')
  write(ifileunit, '(A)')   'structure and fluid field: force_FSI vector'
  write(ifileunit,110) (solid_force_FSI(1,in),solid_force_FSI(2,in),0.0,in=1,nn_solid), &
                       (f_fluids(1,in),f_fluids(2,in),0.0,in=1,nn), &
		       (0.0,0.0,0.0,in=1,nn_inter)
  close(ifileunit)

 !...Write pressure output in ens_movie.pre*
  write(*,*) 'writing... ', name_file2
  open(ifileunit, file=name_file2, form='formatted')
  write(ifileunit, '(A)') 'structure and fluid field: pressure'  
  write(ifileunit,110) (solid_pave(in),in=1,nn_solid),(d(3,in),in=1,nn),(0.0,in=1,nn_inter)
  close(ifileunit)

 !...Write stress output in ens_movie.stress*
  write(*,*) 'writing... ', name_file3
  open(ifileunit, file=name_file3, form='formatted')
  write(ifileunit, '(A)') 'structure field: stress'  
  ! Solid solver tensor order: 1->11,2->22,3->12,4->33
  ! Ensight tensor order: 1->11,2->22,3->33,4->12,5->13,6->23
  write(ifileunit,110) (solid_stress(1,in),solid_stress(2,in),0.0,solid_stress(3,in), &
	                    0.0, 0.0,in=1,nn_solid),  &
                       (fluid_stress(1,in),fluid_stress(2,in),fluid_stress(3,in),fluid_stress(4,in),                 &
                         0.0, 0.0,in=1,nn),&
			(0.0,0.0,0.0,0.0,0.0,0.0,in=1,nn_inter)
  close(ifileunit)

 !...Write strain output in ens_movie.strain*
  write(*,*) 'writing... ', name_file4
  open(ifileunit, file=name_file4, form='formatted')
  write(ifileunit, '(A)') 'structure field: strain'  
  write(ifileunit,110) (solid_strain(1,in),solid_strain(2,in),0.0,solid_strain(3,in),  &
                        0.0, 0.0,in=1,nn_solid),  &
                       (fluid_strain(1,in),fluid_strain(2,in),fluid_strain(3,in),                 &
                        fluid_strain(4,in), 0.0, 0.0,in=1,nn), &
			(0.0,0.0,0.0,0.0,0.0,0.0,in=1,nn_inter)
  close(ifileunit)

  write(*,*) 'writing... ', name_file6
  open(ifileunit, file=name_file6, form='formatted')
  write(ifileunit, '(A)') 'structure and fluid field: I_solid'
  write(ifileunit,110) (solid_pave(in),in=1,nn_solid),(I_solid(in),in=1,nn),(0.5,in=1,nn_inter)
  close(ifileunit)

  write(*,*) 'writing... ', name_file7
  open(ifileunit, file=name_file7, form='formatted')
  write(ifileunit, '(A)') 'structure and fluid field: I_fluid'
  write(ifileunit,110) (solid_pave(in),in=1,nn_solid),(I_fluid(in),in=1,nn),(0.5,in=1,nn_inter)
  close(ifileunit)



endif

  return
end subroutine zfem_ensFluid


end module ensight_output
