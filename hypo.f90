!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   *.fi files are used to shorten hypo.f (keeping the overview)
!   the include command reads these files and replaces the include line
!   with the content of these files

subroutine hypo
  use global_simulation_parameter
  use global_constants
  use run_variables
  use delta_nonuniform
  use solid_variables
  use fluid_variables
  use r_common, only: ninit, vis_solid
  use meshgen_fluid
  use meshgen_solid
  use form
  use ensight_output
  use mpi_variables ! call mpi variable module
  implicit none
  include 'mpif.h'
!==============================	  
! Definition of variables
  integer :: klok,j

  integer infdomain(nn_solid)
  real(8) mass_center(2)
! Variables for fixed solid points on fluid boundary
  integer :: node_sfcon, node_sfcon1  
  integer sfcon_1(125) 
  integer sfcon(250)
  real(8) sfxyz(nsd,250)
  integer inode_sf
! Variables for different fluid density using by implicit form  
  integer mdata(nn_solid)
  integer n_mdata
  real(8) res_l0
  real(8) del_l0

  integer ie, inen
  integer tmp_index(ne)
! For output pressure on the solid nodes
  real(8) pre_inter(nn_solid)
!============================
! Variables for boudary equations
  integer bc4el(ne_inflow) ! 10 is the number of nodes on edge 4
  real(8) res_bc(nsd,nn) ! residual comming from nature B.C. integration 
  real(8) time
  real(8) time_com
!============================
! Define local variables
  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"

!============================
! Define varibales on each processor
  include "hypo_declaration_part.fi"
!===============================================================
! Prepare for calculation, read in inputs or restart information
  include "hypo_restart_file_check.fi"
  include "hypo_prepare_solid.fi"
  include "hypo_prepare_fluid.fi"
!===================================
! Prepare for MPI
!call readpartele(partele)
  include "hypo_prepare_part.fi"
  include "hypo_prepare_com_node.fi"
!=============================
! define the influence domain matrix
 ! integer infdomain(nn_solid)
      call mpi_barrier(mpi_comm_world,ierror)
  write(*,*) 'myid', myid, 'nn_local', nn_local, 'ne_local', ne_local !id for debuger

!=============================
! save the orignal position of solid nodes at fluid boundary
!  do inode_sf=1,node_sfcon
!     sfxyz(1:nsd,inode_sf)=solid_coor_init(1:nsd,sfcon(inode_sf))
!  end do
  vis_solid=vis_liq
if (edge_inflow .ne. 0) then
call edgeele(edge_inflow,rng,neface,ne,bc4el,ne_inflow)
end if

!===================================
! save the orignal position of solid nodes at fluid boundary
  do inode_sf=1,node_sfcon
     sfxyz(1:nsd,inode_sf)=solid_coor_init(1:nsd,sfcon(inode_sf))
  end do

  if (restart == 0) then
	if (myid == 0) then
    	 include 'hypo_write_output.fi'
	end if
  else
     include "hypo_restart_read.fi"
  endif

!=================================================================
!						 time loop	
!=================================================================
  time_loop: do its = nts_start,nts !.....count from 1 or restart-timestep to number of timesteps
      call mpi_barrier(mpi_comm_world,ierror)


	if (myid ==0) then

    	 write (6,*) ' '
    	 write (6,*) 'TIME STEP = ', its
    	 write (6,*) ' '
    	 write (7,*) ' '
    	 write (7,*) 'TIME STEP = ', its
    	 write (7,*) ' '
	
!=================================================================
! Write restart information in binary file

    	 include "hypo_restart_write.fi"
	end if


     tt = tt + dt    !....update real time
     klok = klok + 1 !....update counter for output

	if (myid ==0) then
    	 write (6,'("  physical time = ",f7.3," s")') tt
    	 write (7,'("  physical time = ",f7.3," s")') tt
	end if

! choise of the interpolation method
if (ndelta==1) then

! call the communication and solid solver only on processor 0
	if (myid == 0) then

!=================================================================
! Construction of the dirac deltafunctions at actual solid and fluid node positions
     call delta_initialize(nn_solid,solid_coor_curr,x,ien,dvolume)

!=================================================================
! Solid solver
 write(*,*) 'starting solid solver'
    call solid_solver(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel,  &
                     solid_pave,solid_stress,solid_strain,solid_force_FSI,mtype)

!=================================================================
! Distribution of the solid forces to the fluid domain
!   f^fsi(t)  ->  f(t)
 write(*,*) 'calculating delta'
     call delta_exchange(solid_force_FSI,nn_solid,f_fluids,nn,ndelta,dvolume,nsd,  &
                         delta_exchange_solid_to_fluid)
	endif
!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
      call mpi_barrier(mpi_comm_world,ierror)
      call mpi_bcast(f_fluids(1,1),nsd*nn,mpi_double_precision,0,mpi_comm_world,ierror)
      include "hypo_fluid_solver.fi"


	if (myid == 0) then
!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)
    call delta_exchange(solid_vel,nn_solid,d(1:nsd,:),nn,ndelta,dvolume,nsd, &
					  delta_exchange_fluid_to_solid)
	endif

else if (ndelta==2) then
!=================================================================
! Construction of the FEM influence domain
     call search_inf_pa(solid_coor_curr,x,nn,nn_solid,nsd,ne,nen,ien,infdomain,&
			ne_intlocal,ien_intlocal)

!================================================================
! Merge the auxillary finf array
    call mergefinf(infdomain,nn_solid,mdata,n_mdata)

!=================================================================
! Solid solver
!if (myid ==0) then
! write(*,*) 'starting solid solver'
!end if
    call solid_solver(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel,  &
                     solid_pave,solid_stress,solid_strain,solid_force_FSI,mtype)

!=================================================================
! Set the FSI force for the solid nodes at fluid boundary to be zero
 do inode_sf=1,node_sfcon
   solid_force_FSI(1:nsd,sfcon(inode_sf))=0.0
 end do


!=================================================================
! Distribution of the solid forces to the fluid domain
!   f^fsi(t)  ->  f(t)
if (myid ==0) then
 write(*,*) 'calculating delta'
end if
     call data_exchange_FEM(solid_force_FSI,nn_solid,f_fluids,nn,dvolume,nsd,  &
                         2,ne,nen,ne_solid,nen_solid,&
                        solid_coor_curr,solid_fem_con,x,ien,infdomain,d(nsd+1,:),pre_inter)
!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
      call mpi_barrier(mpi_comm_world,ierror)
time=mpi_wtime()
     include "hypo_fluid_solver.fi"
time=mpi_wtime()-time
if (myid == 0) write(*,*) 'Time for fluid solver', time
!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)
! swith button should be added , right now use 1 or 2 first
    call data_exchange_FEM(solid_vel,nn_solid,d(1:nsd,:),nn,dvolume,nsd, &
			1,ne,nen,ne_solid,nen_solid,&
                       solid_coor_curr,solid_fem_con,x,ien,infdomain,d(nsd+1,:),pre_inter)
end if


!=================================================================
!uPDAte solid domain
!    call solid_update(klok,solid_fem_con,solid_coor_init,solid_coor_curr,  &
!                     solid_vel,solid_prevel,solid_accel)

	include "solid_update_new.fi"



!-----------------------------------------------------------------
! Set the solid nodes at the fluid boundary at their original position
  do inode_sf=1,node_sfcon
     solid_coor_curr(1:nsd,sfcon(inode_sf))=sfxyz(1:nsd,inode_sf)
  end do
!    open(unit=8406, file='masscenter.txt', status='unknown')

!    mass_center(1)=sum(solid_coor_curr(1,:))/nn_solid
!    mass_center(2)=sum(solid_coor_curr(2,:))/nn_solid
!    write(8406,*)  mass_center(:)
!=================================================================
! Volume correction  
!   if (mod(its,10) .eq. 9) then
!   call volcorr(solid_coor_curr,nsd_solid,nen_solid,solid_fem_con,nn_solid,ne_solid, &
!                solid_coor_init)
!   write(*,*) 'volume correction applied'
!   call energy_fluid(x,d(1:nsd,:),ien)
!   end if
!=================================================================
! Write output file every ntsbout steps
!    Out put the fluid pressure at solid nodals
     solid_pave(:)=pre_inter(:)
	if (myid == 0) then
     include "hypo_write_output.fi"
	endif
  enddo time_loop


end subroutine hypo
