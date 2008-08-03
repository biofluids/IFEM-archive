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
  use r_common, only: ninit
  use meshgen_fluid
  use meshgen_solid
  use form
  use ensight_output
  implicit none

!==============================	  
! Definition of variables
  integer :: klok,j
 
  integer infdomain(nn_solid)
 
  real(8) mass_center(2)
  
  integer :: node_sfcon, node_sfcon1
!  integer sfcon_1(219)
!  integer sfcon(438)
!  real(8) sfxyz(nsd,438)
!  integer inode_sf

  integer mdata(nn_solid)
  integer n_mdata

!============================
! Define local variables
  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"
!===============================================================
! Prepare for calculation, read in inputs or restart information

  include "hypo_restart_file_check.fi"
  include "hypo_prepare_solid.fi"
  include "hypo_prepare_fluid.fi"
!=============================
! save the orignal position of solid nodes at fluid boundary
 ! do inode_sf=1,node_sfcon
 !    sfxyz(1:nsd,inode_sf)=solid_coor_init(1:nsd,sfcon(inode_sf))
 ! end do
!-------------------------------------------------------------  
  if (restart == 0) then
     include 'hypo_write_output.fi'
  else
     include "hypo_restart_read.fi"
  endif

!=================================================================
!						 time loop	
!=================================================================
  time_loop: do its = nts_start,nts !.....count from 1 or restart-timestep to number of timesteps

     write (6,*) ' '
     write (6,*) 'TIME STEP = ', its
     write (6,*) ' '
     write (7,*) ' '
     write (7,*) 'TIME STEP = ', its
     write (7,*) ' '

!=================================================================
! Write restart information in binary file

     include "hypo_restart_write.fi"

     tt = tt + dt    !....update real time
     klok = klok + 1 !....update counter for output

     write (6,'("  physical time = ",f9.6," s")') tt
     write (7,'("  physical time = ",f9.6," s")') tt


! choise of the interpolation method
if (ndelta==1) then
!=================================================================
! Construction of the dirac deltafunctions at actual solid and fluid node positions
     call delta_initialize(nn_solid,solid_coor_curr,x,ien,dvolume)

!=================================================================
! Solid solver
 write(*,*) 'starting solid solver'
    call solid_solver(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel,  &
                     solid_pave,solid_stress,solid_strain,solid_force_FSI)

!=================================================================
! Distribution of the solid forces to the fluid domain
!   f^fsi(t)  ->  f(t)
 write(*,*) 'calculating delta'
     call delta_exchange(solid_force_FSI,nn_solid,f_fluids,nn,ndelta,dvolume,nsd,  &
                         delta_exchange_solid_to_fluid)

!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
     include "hypo_fluid_solver.fi"

!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)
    call delta_exchange(solid_vel,nn_solid,d(1:nsd,:),nn,ndelta,dvolume,nsd, &
					  delta_exchange_fluid_to_solid)
else if (ndelta==2) then
!=================================================================
! Construction of the FEM influence domain
     call search_inf(solid_coor_curr,x,nn,nn_solid,nsd,ne,nen,ien,infdomain)

!================================================================
! Merge the auxillary finf array
    call mergefinf(infdomain,nn_solid,mdata,n_mdata)

!=================================================================
! Solid solver
 write(*,*) 'starting solid solver'
    call solid_solver(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel,  &
                     solid_pave,solid_stress,solid_strain,solid_force_FSI)


!=================================================================
! Set the FSI force for the solid nodes at fluid boundary to be zero
! do inode_sf=1,node_sfcon
!   solid_force_FSI(1:nsd,sfcon(inode_sf))=0.0
! end do
!------------------------------------------------------------------




!=================================================================
! Distribution of the solid forces to the fluid domain
!   f^fsi(t)  ->  f(t)
 write(*,*) 'calculating delta', solid_fem_con(1,3)
     call data_exchange_FEM(solid_force_FSI,nn_solid,f_fluids,nn,dvolume,nsd,  &
                         2,ne,nen,ne_solid,nen_solid,&
                        solid_coor_curr,solid_fem_con,x,ien,infdomain)

!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
     include "hypo_fluid_solver.fi"

!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)
! swith button should be added , right now use 1 or 2 first
    call data_exchange_FEM(solid_vel,nn_solid,d(1:nsd,:),nn,dvolume,nsd, &
			1,ne,nen,ne_solid,nen_solid,&
                       solid_coor_curr,solid_fem_con,x,ien,infdomain)
end if

!=================================================================
!uPDAte solid domain
    call solid_update(klok,solid_fem_con,solid_coor_init,solid_coor_curr,  &
                     solid_vel,solid_prevel,solid_accel)


!-----------------------------------------------------------------
! Set the solid nodes at the fluid boundary at their original position
!  do inode_sf=1,node_sfcon
!     solid_coor_curr(1:nsd,sfcon(inode_sf))=sfxyz(1:nsd,inode_sf)
!  end do


    open(unit=8406, file='masscenter.txt', status='unknown')

    mass_center(1)=sum(solid_coor_curr(1,:))/nn_solid
    mass_center(2)=sum(solid_coor_curr(2,:))/nn_solid
    write(8406,*)  mass_center(:)
!===========================================================
! Giving solid coor
!solid_coor_curr(1,1)=0.35d0
!solid_coor_curr(2,1)=0.3d0

!solid_coor_curr(1,2)=0.55d0
!solid_coor_curr(2,2)=0.3d0

!solid_coor_curr(1,3)=0.35d0
!solid_coor_curr(2,3)=0.425d0

!solid_coor_curr(1,4)=0.55d0
!solid_coor_curr(2,4)=0.425d0

!solid_coor_curr(1,5)=0.35d0
!solid_coor_curr(2,5)=0.55d0

!solid_coor_curr(1,6)=0.55d0
!solid_coor_curr(2,6)=0.55d0
!=================================================================
! Volume correction  
!   if (mod(its,100) .eq. 99) then
!   call volcorr(solid_coor_curr,nsd_solid,nen_solid,solid_fem_con,nn_solid,ne_solid, &
!                solid_coor_init)
!   write(*,*) 'volume correction applied'
!   call energy_fluid(x,d(1:nsd,:),ien)
!   end if
!=================================================================
! Write output file every ntsbout steps
     include "hypo_write_output.fi"

  enddo time_loop



end subroutine hypo
