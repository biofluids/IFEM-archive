	subroutine error(string, param, fatal)
c        include "mpif.h"
	character string*(*)
	integer param, strlen
	logical fatal
c        integer myid,comm,ierr

c        call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr)

c	if (myid.eq.0) then
	   if (string.ne.' ') then
	      strlen = len(string)
	      do while (string(strlen:strlen).eq.' ')
		 strlen = strlen - 1
	      end do

	      if (param.ne.-999) then
		 write(unit=0, fmt='(a$)') string
		 write(unit=0, fmt=*) param
	      else
		 write(unit=0, fmt='(a)') string(1:strlen)
	      end if
	   end if
c	end if

	if (fatal) call exit(1)

	return
	end
