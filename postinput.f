c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	parseinput.fcm                                                       c
c	---------------------------------------------------------------------c
c	sets parameters from standard input or control file                  c
c	---------------------------------------------------------------------c
c	910919 - written                                                     c
c	---------------------------------------------------------------------c
c	M. Behr [AHPCRC]                                                     c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine postinput
	  
	  include "global.h"
	  
	  character*32 key, keyaux
	  character*8 onoff
	  logical fctrl, getkey, isatty
	  logical enough
	  data enough /.false./
	  
c  command loop
	  do while (.not.enough)
		enough = .not.getkey(key)
		
		if ((key.eq.'abort').or.(key.eq.'exit')) then
		  call exit(1)
		else if ((key.eq.'done').or.(key.eq.'quit')) then
		  enough = .true.
c  ---------------------------------------------------------------------c
		else if (key.eq.'t_start') then
		  call getreal(t_start)
		else if (key.eq.'dt') then !..........time step
		  call getreal(dt)
		else if (key.eq.'iquad') then
		  call getint(iquad)
		else if (key.eq.'nen') then
		  call getint(nen)
		else if ((key.eq.'ndf').or.(key.eq.'ndof')) then
		  call getint(ndf)
		else if (key.eq.'ne') then
		  call getint(ne)
		else if (key.eq.'nqd') then
		  call getint(nqd)
		else if (key.eq.'maxconn') then
		  call getint(maxconn)
		else if (key.eq.'nn') then
		  call getint(nn)
		else if (key.eq.'nrng') then
		  call getint(nrng)
		else if (key.eq.'nit') then
		  call getint(nit)
		else if (key.eq.'nsd') then !.........space dimensions
		  call getint(nsd)
		else if (key.eq.'nts') then
		  call getint(nts)
		else if (key.eq.'idisk') then
		  call getint(idisk)
		else if (key.eq.'ntsbout') then !.....time steps between disk output
		  call getint(ntsbout)
		else if (key.eq.'surf') then
		  call getint(surf(0))
		  do i=1,surf(0)
			call getint(surf(i))
		  enddo
		else if (key.eq.'hydro') then
		  call getint(hydro)
		else if (key.eq.'restart') then !.....restart
		  call getstr(onoff)
		  restart = (onoff.ne.'off')
		else if (key.eq.'alpha') then !.......time discretization parameter
		  call getreal(alpha)
		else if (key.eq.'delta') then !.......stabilization deltas
		  call getint(idelta)
		  call getreal(delta(idelta))
		else if (key.eq.'iscaling') then !....scaling choice (diag - 1)
		  call getint(iscaling)
		else if (key.eq.'inner') then !....inner GMRES iterations
		  call getint(inner)
		else if (key.eq.'outer') then !....outer GMRES iterations
		  call getint(outer)
		else if (key.eq.'steady') then !......steady state
		  call getstr(onoff)
		  steady = (onoff.ne.'off')
		else if (key.eq.'conserve') then !......steady state
		  call getstr(onoff)
		  conserve = (onoff.ne.'off')
		else if (key.eq.'turb_kappa') then
		  call getreal(turb_kappa)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		else if (key.eq.'vis_liq')then
		  call getreal(vis_liq)
		else if (key.eq.'vis_gas')then
		  call getreal(vis_gas)
		else if (key.eq.'den_liq')then
		  call getreal(den_liq)
		else if (key.eq.'den_gas')then
		  call getreal(den_gas)
		else if (key.eq.'ref_den') then
		  call getreal(ref_den)
		else if (key.eq.'ref_lgt') then
		  call getreal(ref_lgt)
		else if (key.eq.'ref_vel') then
		  call getreal(ref_vel)
		else if (key.eq.'stokes') then !.....restart
		  call getstr(onoff)
		  stokes = (onoff.ne.'off')
		else if (key.eq.'static') then 
		  call getstr(onoff)
		  static = (onoff.ne.'off')
		else if (key.eq.'hg_vol') then 
		  call getstr(onoff)
		  hg_vol = (onoff.ne.'off')
		else if (key.eq.'taudt') then
		  call getstr(onoff)
		  taudt = (onoff.ne.'off')
		  
		else if (key.eq.'gravity') then
		  call getreal(gravity(1))
		  call getreal(gravity(2))
		  call getreal(gravity(3))
		else if (key.eq.'interface') then
		  call getreal(interface(1))
	      call getreal(interface(2))
  	      call getreal(interface(3))
		  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	    else if (key.eq.'slip') then
		  call getint(ibc)
		  do idf=1,ndf
			call getreal(bv(idf,ibc))
			if(abs(bv(idf,ibc)+999.0).gt.1.0e-6) bc(idf,ibc) = 1
		  enddo 
		  call getreal(bvf(ibc))
		  if(abs(bvf(ibc)+999.0).gt.1.0e-6) bcf(ibc) = 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	    else if (key.eq.'fix') then
		  call getint(ibc)
		  do idf=1,ndf
			call getreal(bv(idf,ibc))
			if(abs(bv(idf,ibc)+999.0).gt.1.0e-6) bc(idf,ibc) = 1
		  enddo 
		  call getreal(bvf(ibc))
		  if(abs(bvf(ibc)+999.0).gt.1.0e-6) bcf(ibc) = 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		  
	    else if (key.eq.'initial') then
		  do idf=1,ndf
			call getreal(ic(idf))
		  enddo 
		  call getreal(icf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		else
		  
		  if ((.not.enough).and.(key.ne.' ').and.(key(1:1).ne.'#'))
	1		   call error("parseinput: strange keyword "//key, -999, .false.)
		end if
		
	  end do
	  
c  further defaults
	  if (ntsbout.eq.0) ntsbout = nts + 1
	  if (steady) alpha = 1.0
	  
	  if      (nen.eq.4) then
		etype = tet
		neface = 4
		nnface = 3
	  else if (nen.eq.8) then
		etype = hex
		neface = 6
		nnface = 4
	  end if
	  
	  if (ndf.eq.0) ndf = 4
	  if (nen.eq.0) nen = 8
	  if ( nq.eq.0) nq  = ndf * nn
	  if (nqf.eq.0) nqf = nn
	  twod = .false.
	  if((nn.gt.2*ne).and.(nen.eq.8)) then
		twod = .true.
		hg_vol = .true.
	  endif
	  
c	  nec = (ne - 1) / numproc + 1
c	  maxnec = nec
c	  if (myid.eq.numproc-1) nec = ne - (numproc - 1) * maxnec
c	  nnc = (nn - 1) / numproc + 1
c	  maxnnc = nnc
c	  if (myid.eq.numproc-1) nnc = nn - (numproc - 1) * maxnnc
c	  nqdc = (nqd - 1) / numproc + 1
c	  maxnqdc = nqdc
c	  if (myid.eq.numproc-1) nqdc = nqd - (numproc - 1) * maxnqdc

	return
	end