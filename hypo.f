c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	main.fcm															 c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   *.fi files are used to shorten hypo.f (keeping the overview)
c   the include command reads these files and replaces the include line
c   with the content of these files

	subroutine hypo

	implicit real*8 (a-h,o-z)
cccccccccccccccccccccccccccccc	  
c Definition of variables
	  
	!...FEM fluid solver
	include "global.h" !.......global values for fluid FEM Solver
	include "malloc.h" !.......memory needed to allocate pointers
	
	!...IBM code
	include 'r_common' !.......declares some IBM-variables as common
	include 'main_common' !....same as above 
c	parameter (maxnn_solids=5000)

	real* 8 shrknode(90,mno)
	integer ncnn(mno),cnn(90,mno)
	real*8 totalvel(3),avgvel(3),mom(3),rs(3),xj(3,3),xx(3,9)

	!...defines variables, which are used locally (not common) in hypo.f

	include "declaration_solid.fi"
	include "declaration_fluid.fi"	

cccccccccccccccccccccccccccccc			
c Prepare for calculation, read in inputs 
	include "prepare_solid.fi"
	include "prepare_fluid.fi"

	if (mno .lt. nnd) then
		write(*,*) 'boost maxnn_solids in hypo.f and delta_nonuniform'
		stop
	endif
      include 'write_output.fi'

ccccccccccccccccccccccccccccccc
c time loop

	do its=1,nts !.....count from 1 to number of timesteps
	   istep = its	! time step number
	   iti	 = its

	   tt = tt + dt !......update real time
	   klok = klok + 1 !.... for ibm output
!	   ntsbout=n_step_wr_ib_user_files  ! output steps are the same

	 write (6,*) ' '
	 write (6,*) ' TIME STEP = ', its
	 write (7,*) ' '
	 write (7,*) ' TIME STEP = ', its

cccccccccccccccccccccccccccccccccccccccccc	  
c Construction of the dirac deltafunctions for nonuniform spacing

	 write(*,*) 'calculating RKPM delta function'
	 if (ndelta.eq.1) then	!non-uniform RKPM delta function
	 call delta_nonuniform(shrknode,cnn,ncnn,nnd,coord_pt,xn,
	1	ien,dvolume)
	 endif

cccccccccccccccccccccccccccccccccccccccc
c   2.  Solid solver - membrane element
	 write(*,*) 'solving solids'

	if (nsd_solid .eq.3) then
	 write(*,*) 'solving for 3-d structure'
	 include "solids_solver.fi"
	elseif (nsd_solid .eq. 0) then ! point
	 write(*,*) 'solving for a point'
	 force_pt(1,1)=-xmg1*sdensi
	 force_pt(2,1)=-xmg2*sdensi
	 force_pt(3,1)=-xmg3*sdensi
	endif
cccccccccccccccccccccccccccccccccccccccccccccccccc
c 2.5 Distribution of the solid forces to the fluid domain
c	F(t + dt)  ->  f(t + dt)							 
	 ibuf=2
	 write(*,*) 'distributing the force onto fluids'
	 call diracdelta(force_pt,nnd, coord_pt,
	1	 f_fluids,nn,
     +	ndelta,shrknode,cnn,ncnn,maxconn,dvolume,ibuf)

cccccccccccccccccccccccccccccccccccccccccccccccccc
c 2. FEM Navier-Stokes Solver - calculates v(t+dt),p(t+dt)
	 write(*,*) 'solving fluids'
	 include "fem_fluid_solver.fi"

cccccccccccccccccccccccccccccccccccccccccccccccccc
c 2. Interpolation fluid velocity -> immersed material points
c	   v_fluid(t+dt)  ->  v_solid(t+dt) 	 
		  
	 ibuf=1
	 write(*,*) 'interpolating fluid velocity onto the solids'
	 call diracdelta(vel_pt,nnd, coord_pt,
	1	 d(1:3,:),nn,
     +	ndelta,shrknode,cnn,ncnn,maxconn,dvolume,ibuf)
	 
ccccccccccccccccccccccccccccccccccccccccccccccccc
c	update solid domain
	 write(*,*) 'updating the variables for the structure'
	 include "solids_update.fi"

c	write output file every ntsbout steps

	 if (mod(its,ntsbout).eq.0) then
	    write(*,*) 'generating the output'
       	include 'write_output.fi'
	    if ((n_ibmfem .eq. 0).and.(n_tec_ens .eq. 0)) then
	       write(*, *) 'generate ensight case file'
	       call zibm_enscase(tt, its,ntsbout)
	    endif
	    if((n_ibmfem .eq. 1).and. (n_tec_ens .eq. 0)) then
	       write(*, *) 'generate ensight case file'
	       call zfem_enscase(dt, its,ntsbout)
	    endif
	 endif

	enddo			!....end of time loop 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	write output summary for ensight *.case 
c	n_ibmfem=0 : ibm
c	n_ibmfem=1 : fem
c	n_tec_ens=0: ensight
c	n_tec_ens=1: tecplot
c

	
	!...stops time counting and write output to screen
	naxx2=time()
	write(*,*) naxx1,naxx2,naxx2-naxx1
	
	end 