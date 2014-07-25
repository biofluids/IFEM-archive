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
  use r_common, only: ninit, vis_solid, density_solid, material_type
  use meshgen_fluid
  use meshgen_solid
  use form
  use ensight_output
  use mpi_variables
  use pml_variables
  implicit none
  include 'mpif.h'
!==============================	  
! Definition of variables
  integer :: klok,j

  integer infdomain(nn_solid)
  real(8) mass_center(nsd_solid)
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
! solid 1st type bc temp varaibles
  real(8) omega_tmp
  real(8) theta_tmp

!============================
! Define local variables
  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"

!============================
! Define fluid varibales on each processor
  include "hypo_declaration_part.fi"
! Define solid varibales on each processor
  include "hypo_declaration_part_solid.fi"

!===============================================================
! Prepare for calculation, read in inputs or restart information
  include "hypo_restart_file_check.fi"
  include "hypo_prepare_solid.fi"
  include "hypo_prepare_fluid.fi"
!===================================
! Prepare for MPI Fludi solver
  include "hypo_prepare_part.fi"
  include "hypo_prepare_com_node.fi"

! Prepare for MPISolid solver
  include "hypo_prepare_part_solid.fi"
  include "hypo_prepare_com_node_solid.fi"
!=============================
! define the influence domain matrix
 ! integer infdomain(nn_solid)
      call mpi_barrier(mpi_comm_world,ierror)
!  write(*,*) 'myid', myid, 'nn_local', nn_local, 'ne_local', ne_local !id for debuger
!=============================
  !vis_solid= - vis_liq 
  vis_solid=  1.0
  I_fluid(:)=0.0
  solid_pave(:)=0.0d0
  solid_vel(:,:) = 0.0d0
  solid_accel(:,:) = 0.0d0
  damp_solid = 50.0
  solid_bcvel(:,:) = 0.0
  solid_bcvel_old(:,:) = 0.0
  outedge=2
  fden(:)=0.0
  fvis(:)=0.0
!-----------------------------
! PML auxiliary variable initialization
! DON'T forget putting this in the restart file!!!!!!!
!  qv(:,:)=0.0
!-----------------------------------------------------------------
! calculate \sigma_x and \sigma_y for PML... and, other things lol
call initializePMLparam(x,node_local,nn_local,send_address,ad_length)
include "reallocateVarPML.f90"
!-----------------------------------------------------------------
!=============================

if (edge_inflow .ne. 0) then
    call edgeele(edge_inflow,rng,neface,ne,bc4el,ne_inflow)
endif

!===================================
! save the orignal position of solid nodes at fluid boundary
if (node_sfcon .ne. 0 ) then
    omega_tmp=1.0
    do inode_sf=1,node_sfcon
        sfxyz(1:nsd,inode_sf)=solid_coor_init(1:nsd,sfcon(inode_sf))
!        solid_vel(1,sfcon(inode_sf)) = -(solid_coor_init(2,sfcon(inode_sf)) - 1.0)/0.1*0.1*omega_tmp
!        solid_vel(2,sfcon(inode_sf)) = (solid_coor_init(1,sfcon(inode_sf)) - 4.0)/0.1*0.1*omega_tmp
    enddo
endif

if (restart == 0) then
    if (myid == 0) then
        include 'hypo_write_output.fi'
    endif
else
    include "hypo_restart_read.fi"
endif

!=================================================================
!                          Time Loop
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
        include "hypo_restart_write.fi"
    endif

    tt = tt + dt    !....update real time
    klok = klok + 1 !....update counter for output

    if (myid==0) then
        write (6,'("  physical time = ",f14.10," s")') tt
        write (7,'("  physical time = ",f14.10," s")') tt
    endif

    ! choice of the interpolation method
    ChoiceInterpolation: if (ndelta==1) then
        !--------------------------------
        ! Update the solid coor first
        solid_coor_pre2(:,:) = solid_coor_pre1(:,:)
        solid_coor_pre1(:,:) = solid_coor_curr(:,:)
        solid_stress(:,:) = 0.0d0
        do ie=1,nn_solid
            solid_stress(1:nsd_solid,ie) = solid_pave(ie)
        enddo

        if (node_sfcon .ne. 0 ) then
            do inode_sf=1,node_sfcon
                solid_accel(1,sfcon(inode_sf))= - 2.0*10.0*pi*cos(2.0*pi*10.0*tt)*10.0
                solid_accel(2,sfcon(inode_sf))= 0
            enddo
        endif
        !-------------------------------
        ! correct the curr solid coor by solving the solid mon equations
        call mpi_barrier(mpi_comm_world,ierror)

        StepsSolvingSolid: if (its .gt. 0) then
            if (myid == 0) write(*,*) '=== Fluid solver have converged for', its,'time steps ==='
            time = mpi_wtime()
            id_solidbc(:,:)=1
            OptionSolveSolid: if (material_type .ne. 9) then
                call form_solidid12(id_solidbc,nsd_solid,nn_solid,&
                                    ien_sbc,ne_sbc,nen_solid,ne_solid,solid_fem_con)
                call solve_solid_disp_pa(solid_coor_init,solid_coor_curr,id_solidbc,solid_fem_con,node_sbc, &
                                         solid_coor_pre1,solid_vel,solid_accel,ien_sbc,solid_stress,solid_bcvel,mtype,&
                                         ne_intlocal_solid,ien_intlocal_solid,nn_local_solid,&
                                         node_local_solid,send_address_solid,ad_length_solid,&
                                         global_com_solid,nn_global_com_solid,local_com_solid,nn_local_com_solid)
                if (node_sfcon .ne. 0 ) then
                    do inode_sf=1,node_sfcon
                        solid_coor_curr(1,sfcon(inode_sf))= solid_coor_init(1,sfcon(inode_sf)) - sin(2.0*pi*10.0*tt)*dt 
                        solid_coor_curr(2,sfcon(inode_sf))= solid_coor_init(2,sfcon(inode_sf)) + 0.0d0
                    enddo
                endif
            endif OptionSolveSolid
            call mpi_barrier(mpi_comm_world,ierror)

            time = mpi_wtime()-time
            if (myid == 0)  write(*,*) 'Time for solid solver', time
        endif StepsSolvingSolid

        !-----------------------------------------------------------------
        ! Find the fluid nodes overlapping with solid domain
        if (material_type .ne. 9) then
            time=mpi_wtime()
            call search_inf_re(solid_coor_curr,x,nn,nn_solid,nsd,ne_solid,nen_solid,solid_fem_con,&
                               flag_fnode,node_local,nn_local)
            ! Construction of the dirac deltafunctions at actual solid and fluid node positions
            call rkpm_nodevolume(x,nsd,nn,ien,ne,nen,ien_intlocal,ne_intlocal,dvolume,sp_radius)
            call rkpm_init(solid_coor_curr,nn_solid,x,nsd,nn,dvolume,sp_radius)
            time=mpi_wtime()-time
            if (myid == 0) write(*,*) '---Time for initiate RKPM interpolation function---', time
        endif
        !-----------------------------------------------------------------
        ! Solve Laplace equation get indicatior field
        time=mpi_wtime()
        if (material_type .ne. 9) then
            lp_source(:,:)=0.0
            call source_laplace(x,nn,nsd,solid_coor_curr,nn_solid,solid_fem_con,ne_solid,nen_solid,&
                                ien_sbc,ne_sbc,node_sbc,nn_sbc,lp_source)
            call mpi_barrier(mpi_comm_world,ierror)
            call solve_laplace_pa(lp_source,nsd,nn,nn_solid,ien,ne,nen,x,node_sbc,nn_sbc,I_fluid,flag_fnode,&
                                  ne_intlocal,ien_intlocal,nn_local,node_local,send_address,ad_length,&
                                  global_com,nn_global_com,local_com,nn_local_com)
!            call solve_laplace(lp_source,nsd,nn,nn_solid,ien,ne,nen,x,node_sbc,nn_sbc,I_fluid,flag_fnode)
        endif
        fden(:)=den_liq+density_solid*I_fluid(:)
        fvis(:)=vis_liq+vis_solid*I_fluid(:)
        time=mpi_wtime()-time
        if (myid == 0) write(*,*) '---Time for update indicator field---', time
        !=================================================================
        ! Solid solver
        if ((myid == 0) .and. (material_type .ne. 9)) then
            write(*,*) 'starting my own litte solid solver'
            call res_solid(solid_coor_init,solid_coor_curr,solid_fem_con, solid_force_FSI,&
                           solid_coor_pre1,solid_coor_pre2,id_solidbc,solid_accel,&
                           solid_pave,solid_bcvel,solid_bcvel_old,solid_vel)
!            solid_force_FSI(:,:) =  solid_force_FSI(:,:) + solid_bcforce(:,:)
            ! Set the FSI force for the solid nodes at fluid boundary to be zero
            if (node_sfcon .ne. 0) then
                do inode_sf=1,node_sfcon
                    solid_force_FSI(1:nsd,sfcon(inode_sf))=0.0
                enddo
            endif
            ! Distribution of the solid forces to the fluid domain: f^fsi(t)  ->  f(t)
            write(*,*) 'calculating delta'
            call delta_exchange(solid_force_FSI,nn_solid,f_fluids,nn,ndelta,dvolume,nsd,  &
                                delta_exchange_solid_to_fluid,solid_pave,d(ndf,:))
            do ie = 1, nn
                f_fluids(:,ie) = f_fluids(:,ie) * fden(ie)
            enddo
        elseif (material_type == 9) then
            f_fluids(:,:) = 0.0
        endif
        !=================================================================
        call mpi_barrier(mpi_comm_world,ierror)
        call mpi_bcast(f_fluids(1,1),nsd*nn,mpi_double_precision,0,mpi_comm_world,ierror)
!        call mpi_bcast(fden(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
!        call mpi_bcast(fvis(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
!        call mpi_bcast(I_fluid(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
        time=mpi_wtime()
!------------------------------------------------------------------------
!if (tt==0.0) then
!    steady = .true.
!    include "hypo_fluid_solver.f90"
!    steady = .false.
!endif
!------------------------------------------------------------------------
        include "hypo_fluid_solver.f90"

        time=mpi_wtime()-time
        if (myid == 0) write(*,*) '---Time for fluid solver---', time
        !=================================================================
        ! Interpolation fluid velocity -> immersed material points: v^f(t+dt)  ->  v^s(t+dt)
        solid_bcvel_old(:,:) =  solid_bcvel(:,:)
        if (material_type .ne. 9) then
            call delta_exchange(solid_bcvel,nn_solid,d(1:nsd,:),nn,ndelta,dvolume,nsd, &
                                delta_exchange_fluid_to_solid,solid_pave,d(ndf,:))
        endif

    elseif (ndelta==2) then
        !=================================================================
        ! Not an option now
        if (myid == 0) write(*,*) 'Not an option now as ndelta == 2 for semi-implicit FSI'
        stop
    endif ChoiceInterpolation 
    if (myid == 0) then
        include "hypo_write_output.fi"
    endif
enddo time_loop

end subroutine hypo