      subroutine timer(p,oldp)
      include "global.h"

      integer p,oldp

c parameters for the time
c 1.  setup time
c 2.  communication time
c 3.  input & output
c 4.  boundary setup
c 5.  gmres
c 6.  total time

      if (totaltime(p).eq.0) then
         starttime(p) = MPI_WTIME() 
      endif
      if (p.eq.oldp) then
         endtime(p) = MPI_WTIME()
         totaltime(p) = totaltime(p) + (endtime(p) - starttime(p))
      else 
         starttime(p) = MPI_WTIME()
      endif
      oldp = p
      end
