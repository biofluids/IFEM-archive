c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readx(xn)

	include "global.h"
	real* 8 xn(nsd,nnc)
      integer lock,ierr,status(MPI_STATUS_SIZE)

	maxrecl = nsd*maxnnc*8
      lock = 1
	if (myid.gt.0) then
      call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
     &                MPI_COMM_WORLD,status,ierr)
	end if
	call ewd_open("mxyz", ifp)
	call ewd_lseek(ifp, myid*maxrecl,0)
	call ewd_read(ifp, xn, nsd*nnc*8)
	call ewd_close(ifp)
	if (myid.lt.numproc-1) then
      call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
     &                MPI_COMM_WORLD,ierr)
	end if

      return
      end
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readien(ien)

	include "global.h"
	integer ien(nen,nec), ienieee(nen,nec/2+1)
      integer lock,ierr,status(MPI_STATUS_SIZE)

	maxrecl = nen * maxnec * 4
      lock = 1
	if (myid.gt.0) then
      call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
     &                MPI_COMM_WORLD,status,ierr)
	end if
	call ewd_open("mien", ifp)
	call ewd_lseek(ifp, myid*maxrecl, 0)
	call ewd_read(ifp, ienieee, nen*nec*4)
	call ewd_close(ifp)
	if (myid.lt.numproc-1) then
      call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
     &                MPI_COMM_WORLD,ierr)
	end if

      call iei2cray(ienieee,ien,nen*nec)

	return
	end
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readrng(rngface)

	include "global.h"
	integer rngface(neface,nec), rngfaceieee(neface,nec/2+1)
      integer lock,ierr,io,status(MPI_STATUS_SIZE)

	maxrecl = neface * maxnec * 4
      lock = 1
	if (myid.gt.0) then
      call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
     &                MPI_COMM_WORLD,status,ierr)
	end if
	call ewd_open("mrng", ifp)
	call ewd_lseek(ifp, myid*maxrecl, 0)
	call ewd_read(ifp, rngfaceieee, neface*nec*4)
	call ewd_close(ifp)
      if (myid.lt.numproc-1) then
      call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
     &                MPI_COMM_WORLD,ierr)
	end if

      call iei2cray(rngfaceieee,rngface,neface*nec)

	do ieface=1,neface
	do iec=1,nec
	  if(rngface(ieface,iec).lt.0) rngface(ieface,iec) = 0
	enddo
	enddo

	mynrng = 0
	do ieface=1,neface
	do iec=1,nec
		mynrng = max(mynrng, rngface(ieface,iec))
	end do
	end do

	if (myid.eq.0) then
		do io=1,numproc-1
                call MPI_RECV(nrng,1,MPI_INTEGER,io,101,
     &                        MPI_COMM_WORLD,status,ierr)
		mynrng = max(mynrng, nrng)
		end do
		nrng = mynrng
	      write(6,1000) nrng
 1000		format("mrng: ",i2," boundaries")
	else
            call MPI_SEND(mynrng,1,MPI_INTEGER,0,101,
     &                        MPI_COMM_WORLD,ierr)

	end if

        call MPI_BCAST(nrng,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	return
	end