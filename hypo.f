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
	real* 8 newcoor(nsd,nn+801) 
	parameter (maxnn_solids=1000)

	real* 8 shrknode(90,maxnn_solids)
	integer ncnn(maxnn_solids),cnn(90,maxnn_solids)

	!...defines variables, which are used locally (not common) in hypo.f

	include "declaration_solid.fi"
	include "declaration_fluid.fi"	

cccccccccccccccccccccccccccccc			
c Prepare for calculation, read in inputs 
	include "prepare_solid.fi"
	include "prepare_fluid.fi"
	if (maxnn_solids .lt. nnd) then
		write(*,*) 'boost maxnn_solids in hypo.f and delta_nonuniform'
		stop
	endif

	newcoor(1,:)=coor(:,1)
	newcoor(3,:)=coor(:,2)
	newcoor(2,:)=0.0d0
	newcoor(:,nnd+1:nnd+nn)=xn(:,:)

ccccccccccccccccccccccccccccccc
c time loop

	do its=1,nts !.....count from 1 to number of timesteps
	   istep = its	! time step number
	   iti	 = its

	   tt = tt + dt !......update real time
	   klok = klok + 1 !.... for ibm output
	   ntsout=n_step_wr_ib_user_files  ! output steps are the same

	 write (6,*) ' '
	 write (6,*) ' TIME STEP = ', its
	 write (7,*) ' '
	 write (7,*) ' TIME STEP = ', its

cccccccccccccccccccccccccccccccccccccccccc	  
c Construction of the dirac deltafunctions for nonuniform spacing

	 write(*,*) 'calculating RKPM delta function'
	 if (ndelta.eq.1) then	!non-uniform RKPM delta function
	 call delta_nonuniform(shrknode,cnn,ncnn,nnd,coord_pt,xn,
	1	ien)
	 endif

cccccccccccccccccccccccccccccccccccccccc
c   2.  Solid solver - membrane element
	 write(*,*) 'solving solids'
	 include "solids_solver.fi"
	   
cccccccccccccccccccccccccccccccccccccccccccccccccc
c 2.5 Distribution of the solid forces to the fluid domain
c	F(t + dt)  ->  f(t + dt)							 
	 ibuf=2
c	 write(*,*) 'distributing the force onto fluids'
c	 call diracdelta(force_pt,nnd, coord_pt,
c	1	 f_fluids,nn,
c     +	ndelta,shrknode,cnn,ncnn,maxconn,ibuf)

cccccccccccccccccccccccccccccccccccccccccccccccccc
c 2. FEM Navier-Stokes Solver - calculates v(t+dt),p(t+dt)
	 write(*,*) 'solving fluids'
	 include "fem_fluid_solver.fi"

cccccccccccccccccccccccccccccccccccccccccccccccccc
c 2. Interpolation fluid velocity -> immersed material points
c	   v_fluid(t+dt)  ->  v_solid(t+dt) 	 
		  
	 ibuf=1
	 write(*,*) 'interpolating fluid velocity onto the solids'
c	 call diracdelta(vel_pt,nnd, coord_pt,
c	1	 d(1:3,:),nn,
c     +	ndelta,shrknode,cnn,ncnn,maxconn,ibuf)
	 
ccccccccccccccccccccccccccccccccccccccccccccccccc
c	update solid domain
	 write(*,*) 'updating the variables for the structure'
	 include "solids_update.fi"

c	write output file every ntsout steps

	 if (mod(its,ntsout).eq.0) then
	    write(*,*) 'generating the output'
       	include 'write_output.fi'
	 endif
	
	newcoor(1,:)=coord_pt(1,:)
	newcoor(2,:)=coord_pt(2,:)
	newcoor(3,:)=coord_pt(3,:)
	newcoor(:,nnd+1:nnd+nn)=xn(:,:)
	enddo			!....end of time loop 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	write output summary for ensight *.case 
c	n_ibmfem=0 : ibm
c	n_ibmfem=1 : fem
c	n_tec_ens=0: ensight
c	n_tec_ens=1: tecplot
c
	if ((n_ibmfem .eq. 0).and.(n_tec_ens .eq. 0)) then
	   write(*, *) 'generate ensight case file'
	   call zibm_enscase(time_value, currentstep)
	end if
	if((n_ibmfem .eq. 1).and. (n_tec_ens .eq. 0)) then
	   write(*, *) 'generate ensight case file'
	   call zfem_enscase(time_value, currentstep)
	end if
	
	!...stops time counting and write output to screen
	naxx2=time()
	write(*,*) naxx1,naxx2,naxx2-naxx1
	
	end 