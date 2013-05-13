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
  use r_common, only: ninit, vis_solid, density_solid
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
  real(8) sfxyz(nsd,node_sfcon)
  integer inode_sf
! Variables for different fluid density using by implicit form  
  integer mdata(nn_solid)
  integer n_mdata
  real(8) res_l0
  real(8) del_l0

  integer ie, inen
! For output pressure on the solid nodes
  real(8) pre_inter(nn_solid)
!============================
! Variables for boudary equations
  integer bc4el(ne_inflow) ! 10 is the number of nodes on edge 4
  real(8) res_bc(nsd,nn) ! residual comming from nature B.C. integration 
  real(8) time
  real(8) time_com
  real(8) pin_s
!============================
character(len=14) filename
character(len=14) resfilename
character(len=7) fileroot
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
!  vis_solid= - vis_liq 
  vis_solid=  1.0
  I_fluid(:)=0.0
  solid_pave(:)=0.0d0
  damp_solid = 10.0
  solid_vel(:,:) = 0.0
  solid_bcvel(:,:) = 0.0
  solid_bcvel_old(:,:) = 0.0
!=============================

if (edge_inflow .ne. 0) then
call edgeele(edge_inflow,rng,neface,ne,bc4el,ne_inflow)
end if

!===================================
! save the orignal position of solid nodes at fluid boundary
if (node_sfcon .ne. 0 ) then
  do inode_sf=1,node_sfcon
     sfxyz(1:nsd,inode_sf)=solid_coor_init(1:nsd,sfcon(inode_sf))
  end do
end if

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
    	 write (6,'("  physical time = ",f14.10," s")') tt
    	 write (7,'("  physical time = ",f14.10," s")') tt
	end if
! choise of the interpolation method
if (ndelta==1) then
!--------------------------------
! Update the solid coor first
        solid_coor_pre2(:,:) = solid_coor_pre1(:,:)
        solid_coor_pre1(:,:) = solid_coor_curr(:,:)
!        solid_prevel(1:nsd_solid,1:nn_solid) = solid_vel(1:nsd_solid,1:nn_solid)
	solid_stress(:,:) = 0.0d0
	do ie=1,nn_solid
		solid_stress(1:nsd_solid,ie) = solid_pave(ie)
	end do

!-------------------------------
! correct the curr solid coor by solving the solid mon equations
call form_solidid12(id_solidbc,nsd_solid,nn_solid,ien_sbc,ne_sbc,nen_solid,ne_solid,solid_fem_con)
call solve_solid_disp(solid_coor_init,solid_coor_curr,id_solidbc,solid_fem_con,node_sbc, &
                        solid_coor_pre1,solid_vel,solid_accel,ien_sbc,solid_stress,solid_bcvel,mtype)

!--------------------------------
! Find the fluid nodes overlapping with solid domain
call search_inf_re(solid_coor_curr,x,nn,nn_solid,nsd,ne_solid,nen_solid,solid_fem_con,&
                flag_fnode,node_local,nn_local)

! Construction of the dirac deltafunctions at actual solid and fluid node positions
time=mpi_wtime()
call rkpm_nodevolume(x,nsd,nn,ien,ne,nen,ien_intlocal,ne_intlocal,dvolume,sp_radius)
call rkpm_init(solid_coor_curr,nn_solid,x,nsd,nn,dvolume,sp_radius)

time=mpi_wtime()-time
if (myid == 0) write(*,*) '---Time for initiate RKPM interpolation function---', time
!==================================================================
! Solve Laplace equation get indicatior field
if (myid == 0) then
	time=mpi_wtime()
	lp_source(:,:)=0.0
        call source_laplace(x,nn,nsd,solid_coor_curr,nn_solid,solid_fem_con,ne_solid,nen_solid,&
			ien_sbc,ne_sbc,node_sbc,nn_sbc,lp_source)
	call solve_laplace(lp_source,nsd,nn,nn_solid,ien,ne,nen,x,node_sbc,nn_sbc,I_fluid,flag_fnode)
! Update density and viscosity field for fluid 
!I_fluid(:)=0.0
	fden(:)=den_liq+density_solid*I_fluid(:)
	fvis(:)=vis_liq+vis_solid*I_fluid(:)
	time=mpi_wtime()-time
	if (myid == 0) write(*,*) '---Time for update indicator field---', time
!=================================================================
! Solid solver

	 write(*,*) 'starting my own litte solid solver'
call res_solid(solid_coor_init,solid_coor_curr,solid_fem_con, solid_force_FSI,&
              solid_coor_pre1,solid_coor_pre2,id_solidbc,solid_accel,solid_pave,solid_bcvel,solid_bcvel_old,solid_vel)
!=================================================================
! Set the FSI force for the solid nodes at fluid boundary to be zero
	if (node_sfcon .ne. 0) then
	 do inode_sf=1,node_sfcon
	   solid_force_FSI(1:nsd,sfcon(inode_sf))=0.0
	 end do
	end if 
!=================================================================
! Distribution of the solid forces to the fluid domain
!   f^fsi(t)  ->  f(t)
	 write(*,*) 'calculating delta'
	     call delta_exchange(solid_force_FSI,nn_solid,f_fluids,nn,ndelta,dvolume,nsd,  &
	                         delta_exchange_solid_to_fluid,solid_pave,d(ndf,:))

do ie = 1, nn
f_fluids(:,ie) = f_fluids(:,ie) * fden(ie)
end do
!f_fluids=0.0

endif
!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
      call mpi_barrier(mpi_comm_world,ierror)
      call mpi_bcast(f_fluids(1,1),nsd*nn,mpi_double_precision,0,mpi_comm_world,ierror)
      call mpi_bcast(fden(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
      call mpi_bcast(fvis(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
      call mpi_bcast(solid_pave(1),nn_solid,mpi_double_precision,0,mpi_comm_world,ierror)
      call mpi_bcast(I_fluid(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)

time=mpi_wtime()

      include "hypo_fluid_solver.fi"

time=mpi_wtime()-time
if (myid == 0) write(*,*) '---Time for fluid solver---', time



!	if (myid == 0) then
!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)


!write(*,*) '***solid velocity interpolation is commended out***'

solid_bcvel_old(:,:) =  solid_bcvel(:,:)

    call delta_exchange(solid_bcvel,nn_solid,d(1:nsd,:),nn,ndelta,dvolume,nsd, &
					  delta_exchange_fluid_to_solid,solid_pave,d(ndf,:))


!	endif

else if (ndelta==2) then
!=================================================================
! Not an option now
if (myid == 0) write(*,*) 'Not an option now as ndelta == 2 for semi-implicit FSI'

stop

end if


!=================================================================
!uPDAte solid domain
! Already updated at the beginning of the time step
	
!    call solid_update(klok,solid_fem_con,solid_coor_init,solid_coor_curr,  &
!                     solid_vel,solid_prevel,solid_accel)


!	include "solid_update_new.fi"


!----------------------------------------
!---------------------------------------
!end if ! debug solid disp solver only
!----------------------------------------
!---------------------------------------

!=================================================================
! Write output file every ntsbout steps
!    Out put the fluid pressure at solid nodals
if (node_sfcon .ne. 0) then
  do inode_sf=1,node_sfcon
     solid_coor_curr(1:nsd,sfcon(inode_sf))=sfxyz(1:nsd,inode_sf)
  end do
end if
	if (myid == 0) then
     include "hypo_write_output.fi"
	endif
  enddo time_loop


end subroutine hypo
