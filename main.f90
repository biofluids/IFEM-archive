!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	main.fcm                                                             c
!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   initialization for the fluid and general variables is done in main.f
!   initialization for the solid variables (and input file reading) in hypo.f
program main
  implicit none
      
 !...set standard values 
  call initialize		! for fluids

 !...read configuration files	  
  call parseinput_fluid
  call nondimension
 !...echos input data
  call echoinput

 !...switch to main routine	   
  call hypo

end program main