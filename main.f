c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	main.fcm                                                             c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program main
         include "global.h"

         integer ierr
	   real starttime,endtime,totaltime

            starttime=TIME()
            print *,' starttime ',starttime
	
            call initialize
	    call parseinput
	    call nondimension
	    call echoinput

         call hypo
   
         endtime=TIME()
         totaltime=endtime-starttime
         print *,' endtime ',endtime
         print *,' total   ',totaltime
	  
         end
	  

















