c	ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	parseinput.fcm                                                      c
c	--------------------------------------------------------------------c
c	sets parameters from standard input or control file                 c
c	--------------------------------------------------------------------c
c	general parameter, program control                                  c
c	ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine parseinput
	implicit none

	include "global.h"
	
	character*32 key, keyaux,test
	character*8 onoff
	logical fctrl, getkey, isatty
	logical enough
	data enough /.false./

c Lucy added idf in the local variable as a counter
	integer idf

c       command loop

	do while (.not.enough)
	   enough = .not.getkey(key)
	   if ((key.eq.'abort').or.(key.eq.'exit')) then
	      call exit(1)
	   else if ((key.eq.'done').or.(key.eq.'quit')) then
	      enough = .true.
c---------------------------------------------------------------------c
c---- time loop control
	   else if (key.eq.'t_start') then !.........physical timestep size
	      call getreal(t_start)
	   else if (key.eq.'dt') then !..............time step
	      call getreal(dt)
	   else if (key.eq.'nit') then !.............number of iterations
	      call getint(nit)
	   else if (key.eq.'nsd') then !.............space dimensions
	      call getint(nsd)
	   else if (key.eq.'nts') then !.............number of time steps
	      call getint(nts)
	   else if (key.eq.'restart') then !.........restart
	      call getstr(onoff)
	      restart = (onoff.ne.'off')
	   else if (key.eq.'steady') then !..........steady state
	      call getstr(onoff)
	      steady = (onoff.ne.'off')
	   else if (key.eq.'conserve') then !........steady state
	      call getstr(onoff)
	      conserve = (onoff.ne.'off')
	      
c---- output control
	   else if (key.eq.'idisk') then !...........???
	      call getint(idisk)
	   else if (key.eq.'ntsbout') then !.........time steps between disk output
	      call getint(ntsbout)
	   else if (key.eq.'surf') then !............???
c	      call getint(surf(0))
c	      do i=1,surf(0)
c		     call getint(surf(i))
c	      end do
c	   else if (key.eq.'ntsoutput_ibm') then !...time steps between disk output (old ibm-code)
c	      call getint(n_step_wr_ib_user_files)

c---- physical parameters	      
	   else if (key.eq.'gravity') then !.........acceleration due to gravity
	      call getreal(gravity(1))
	      call getreal(gravity(2))
	      call getreal(gravity(3))
	   else if (key.eq.'interface') then !.......???
	      call getreal(interface(1))
	      call getreal(interface(2))
  	      call getreal(interface(3))
	      
	   else if (key.eq.'initial') then !.........???
	      do idf=1,ndf
		 call getreal(ic(idf))
	      enddo 
	      call getreal(icf)
	   elseif (key.eq.'maxconn') then ! define maximum of influence pts
	      call getint(maxconn)
	      	      
c old ibm stuff	      
c	   else if (key.eq.'ibmfem') then !...........switch between solid model: fem(1)  ibm(0)
c	      call getint(n_ibmfem)	  
c	   else if (key.eq.'tec_ens') then !..........switch between tecplot(1) and ensight output(0)
c	      call getint(n_tec_ens)	  	      
c	   else if (key.eq.'dispforce') then !........switch on(1) or off(0) disturbance
c	      call getint(n_dispforce)	 	      
c	   else if (key.eq.'matlab') then !...........switch on(1) or off(0) matlab output
c	      call getint(n_matlab)	 	      
	      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	   else
	      
	      if ((.not.enough).and.(key.ne.' ').and.(key(1:1).ne.'#'))
	1	call error("parseinput: strange keyword "//key, -999, .false.)
	   end if
	   
	end do
	
c---- further defaults
	if (ntsbout.eq.0) ntsbout = nts + 1

c---- setting equivalent ibm-code values
c	xssf    = 0.5d-6	!.......???
c	xfactor = 1.0d0		!.......???
c	ntf     = 9988		!.......???
c	n_step_run=nts		!.......number of time steps
c	n_step_wr_pc_file=n_step_run

	return
	end
