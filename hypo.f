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

	 write (6,*) ' '
	 write (6,*) ' TIME STEP = ', its
	 write (7,*) ' '
	 write (7,*) ' TIME STEP = ', its

cccccccccccccccccccccccccccccccccccccccc
c   2.  Solid solver - membrane element
	 write(*,*) 'solving solids'

	 include "solids_solver.fi"

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
	
	!...stops time counting and write output to screen
	naxx2=time()
	write(*,*) naxx1,naxx2,naxx2-naxx1
	
	end 