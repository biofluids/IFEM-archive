	subroutine initialize
	use run_variables
      use fluid_variables
	implicit none

      integer :: isd,i,j,idel

	!real*8 :: mass

	call facemap

	alpha = 1.0
	tt = 0.0
	dt = 1.0
	t_start = 0.0
	ref_lgt = 1.0
	ref_vel = 1.0
	ref_den = 1.0
	vis_liq = 1.0
	den_liq = 1.0
	vis_gas = 1.0
	den_gas = 1.0

	do isd = 1,nsd
	   gravity(isd) = 0.0
	   interface(isd) = -999.0
	enddo
	
	turb_kappa = 0.41

	hydro = 0 
	do i=0,21
	   surf(i) = 0
	enddo

	ndf = 4
	nsd = 3
	nn = 0
	ne = 0
	nq = 0
	nqf = 0
	nec = 0
	nrng = 0

	iquad = 2
	nts = 1
	nit = 1
	ntsbout = 1
	idisk = 0

	nfsurf = 0
	fsurf(:) = 0
	calcforce = .false.

	inner = 1
	outer = 1
	iscaling = 1

	restart  = .false.
	stokes  = .false.
	steady  = .false.
	hg_vol  = .false.
	static  = .false.
	taudt  = .false.
	conserve = .false.

	do j=1,21
	   do i=1,ndfpad
	      bc(i,j) = 0
	      bv(i,j) = -999.0
	   enddo
	   bcf(j) = 0
	   bvf(j) = -999.0
	enddo

	do i=1,ndfpad
	   ic(i) = -999.0
	enddo
	icf = -999

	do idel=0,21
	   delta(idel) = 0.0
	end do

	!mass = 1.0

	return
	end


