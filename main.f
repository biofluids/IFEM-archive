c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	main.fcm                                                             c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   initialization for the fluid and general variables is done in main.f
c   initialization for the solid variables (and input file reading) in hypo.f
	program main

	implicit real*8 (a-h,o-z)
c    include "global.h"
	include "r_common"
	include "main_common"
      
c    set standard values 
	call initialize		! for fluids

c    read configuration files	  
	call parseinput_fluid
	call nondimension
c    echos input data
	call echoinput

c    switch to main routine	   
	call hypo
	stop
	end