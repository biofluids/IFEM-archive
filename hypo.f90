!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	main.fcm															 c
!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   *.fi files are used to shorten hypo.f (keeping the overview)
!   the include command reads these files and replaces the include line
!   with the content of these files

subroutine hypo
  use global_simulation_parameter
  use delta_nonuniform
  use solid_variables
  use fluid_variables
  use r_common

  use form
  implicit none

!ccccccccccccccccccccccccccccc	  
! Definition of variables
	  
 !...FEM fluid solver
  include "malloc.fi" !.......memory needed to allocate pointers

  integer :: i,j,k,ie,i_solid,lx,ly,lz
  integer :: nos,ni,noj,klok
  integer :: ntem,node !,num_fiber
  integer :: naxx1,naxx2

  integer time
  external time

  real*8 :: todet,wp,dum,viter,tot_vol,tot_force


  real*8 totalvel(3),avgvel(3),mom(3),rs(3),xj(3,3),xx(3,9)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define local variables

  !include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
! Prepare for calculation, read in inputs 
  include "hypo_prepare_solid.fi"
  include "hypo_prepare_fluid.fi"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write output for initial configuration
  include 'hypo_write_output.fi'




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! time loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  time_loop: do its=1,nts !.....count from 1 to number of timesteps

	 tt = tt + dt    !....update real time
	 klok = klok + 1 !....update counter for output

	 write (6,*) ' '
	 write (6,*) ' TIME STEP = ', its
	 write (6,'(" physical time = ",f7.3," s")') tt
	 write (7,*) ' '
	 write (7,*) ' TIME STEP = ', its
	 write (7,'(" physical time = ",f7.3," s")') tt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construction of the dirac deltafunctions for nonuniform spacing

	 if (ndelta.eq.1) then	!non-uniform RKPM delta function
	    call delta_initialize(nn_solid,solid_coor_curr,xn,ien,dvolume)
	 endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   2.  Solid solver - membrane element

	 write(*,*) 'solving solids'
	 select case (nsd_solid)
	 case (3) ! 3D structure
	    include "hypo_solid_solver.fi"
	 case (0) ! 0D point
	    write(*,*) 'solving for a point'
	    solid_force_FSI(1,1)=-xmg1*density_solid
	    solid_force_FSI(2,1)=-xmg2*density_solid
	    solid_force_FSI(3,1)=-xmg3*density_solid
	 end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Distribution of the solid forces to the fluid domain
!	f^fsi(t + dt)  ->  f(t + dt)
			
	 write(*,*) 'distributing the force onto fluids'
	 call delta_exchange(solid_force_FSI,nn_solid,f_fluids,nn,ndelta,dvolume,delta_exchange_solid_to_fluid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)

	 include "hypo_fluid_solver.fi"
	 write(*,'(" maximum fluid velocity (x dir) = ",f8.5)'),maxval(d(1,:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolation fluid velocity -> immersed material points
!	   v^f(t+dt)  ->  v^s(t+dt) 	 

	 call delta_exchange(solid_vel,nn_solid,d(1:3,:),nn,ndelta,dvolume,delta_exchange_fluid_to_solid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update solid domain
	 
	 include "hypo_solid_update.fi"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write output file every ntsbout steps

	 if (mod(its,ntsbout).eq.0) then
	    write(*,*) "generating the output"
       	include "hypo_write_output.fi"
	    if(n_tec_ens == 0) then
	       call zfem_ensCase(dt, its,ntsbout)
	    endif
	 endif

  enddo	time_loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	write output summary for ensight *.case 
!	n_ibmfem=0 : ibm
!	n_ibmfem=1 : fem
!	n_tec_ens=0: ensight
!	n_tec_ens=1: tecplot
!

	
 !...stops time counting and write output to screen
  naxx2=time()
  write(*,*) naxx1,naxx2,naxx2-naxx1
	
end subroutine hypo