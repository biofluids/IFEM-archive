!========================
!	hypo.f90									
!========================
!   *.fi files are used to shorten hypo.f (keeping the overview)
!   the include command reads these files and replaces the include line
!   with the content of these files

!=======================================
subroutine hypo


!ccccccccccccccccccccccccccccc	  
! Definition of variables
  use global_simulation_parameter
  use run_variables
  use delta_nonuniform
  use solid_variables
  use fluid_variables
  use r_common, only: xg,wgt,xg_tetra,wgt_tetra,du
  use form
  use ensight_output
  implicit none
  !include "malloc.fi" !.......memory needed to allocate pointers FEM fluid solver

  integer :: klok

  integer :: time_start,time_stop

  integer time
  external time



!===================================
! Define local variables

  !include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"	

!============================================			
! Prepare for calculation, read in inputs 
  include "hypo_prepare_solid.fi"
  include "hypo_prepare_fluid.fi"

!============================================
! Write output for initial configuration
  its = 0
  klok = nrestart
  include 'hypo_write_output.fi'


  time_start=time()

!================================================================
!                time loop  
!================================================================
  time_loop: do its=nrestart+1,nrestart+nts !.....count from 1 to number of timesteps

	 tt = tt + dt    !....update real time
	 klok = klok + 1 !....update counter for output

	 write (6,*) ' '
	 write (6,*) 'TIME STEP = ', its
	 write (6,'("  physical time = ",f7.3," s")') tt
	 write (7,*) ' '
	 write (7,*) 'TIME STEP = ', its
	 write (7,'("  physical time = ",f7.3," s")') tt

!==================================================================
! Construction of the dirac deltafunctions for nonuniform spacing

     call delta_initialize(nn_solid,solid_coor_curr,xn,ien,dvolume)

!==================================================================
!   2.  Solid solver

     call solid_solver

!======================================================
! Distribution of the solid forces to the fluid domain
!	f^fsi(t + dt)  ->  f(t + dt)
			
	 call delta_exchange(solid_force_FSI,nn_solid,f_fluids,nn,ndelta,dvolume,delta_exchange_solid_to_fluid)

!===============================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)

	 include "hypo_fluid_solver.fi"
	 
!===========================================================
! Interpolation fluid velocity -> immersed material points
!	   v^f(t+dt)  ->  v^s(t+dt) 	 

     call delta_exchange(solid_vel,nn_solid,d(1:3,:),nn,ndelta,dvolume,delta_exchange_fluid_to_solid)

!==================================
! Update solid domain
	 
	 call solid_update(klok)

!=========================================
! Write output file every ntsbout steps

     include "hypo_write_output.fi"


  enddo	time_loop

	
!...stops time counting and write output to screen
  time_stop=time()
  write(*,*) time_start,time_stop,time_stop-time_start
	
end subroutine hypo