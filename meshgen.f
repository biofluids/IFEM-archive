c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readx(xn)

	include "global.h"
	real* 8 xn(nsd,nn)
c	integer lock,ierr,status(MPI_STATUS_SIZE)
	integer file,offset,endset,i
c	real* 8 x_temp(nsd,nn),xn_temp(nsd,nn)

	file=23
	open(file, FILE="mxyz.dat", STATUS="old")

	do i=1,nn
	   read(file,*) xn(1:3,i)
	enddo	
	close(file)
	return
	end
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readien(ien)
	
	include "global.h"
	integer ien(nen,ne)
c	integer lock,ierr,status(MPI_STATUS_SIZE)
	character*4 ifp
	integer file,offset,endset,ien_temp(nen,ne)

	file=21
	open(file, FILE="mien.dat", STATUS="old")
	do i=1,ne
	   read(file,*) ien(:,i)
	enddo
	close(file)
	return
	end
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readrng(rngface)

	include "global.h"
	integer rngface(neface,ne), rngfaceieee(neface,ne/2+1)
c	integer lock,ierr,io,status(MPI_STATUS_SIZE)
	character*4 ifp
	integer file,offset,endset,rng_temp(neface,ne)

	file=26
	open(file, FILE="mrng.dat", STATUS="old")
	
	do i=1,ne
	   read(file,*) rngface(:,i)
	enddo

	do ieface=1,neface
	   do iec=1,ne
	      if(rngface(ieface,iec).lt.0) rngface(ieface,iec) = 0
	   enddo
	enddo
	
	mynrng = 0
	do ieface=1,neface
	   do iec=1,ne
	      mynrng = max(mynrng, rngface(ieface,iec))
	   end do
	end do
	nrng=mynrng
	close(file)
	return
	end
