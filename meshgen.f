c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readx(xn)

	include "global.h"
	real* 8 xn(nsd,nnc)
	integer lock,ierr,status(MPI_STATUS_SIZE)
	integer file,offset,endset,i
	real* 8 x_temp(nsd,nn),xn_temp(nsd,nnc)

	file=23
	open(file, FILE="mxyz.dat", STATUS="old")
	do i=1,nn
	   read(file,*) x_temp(:,i)
	enddo	

	lock = 1
	if (myid.gt.0) then
	   call MPI_RECV(lock,1,MPI_DOUBLE_PRECISION,myid-1,101,
	1	MPI_COMM_WORLD,status,ierr)
	end if

	offset = myid * maxnnc + 1
	endset = offset + nnc-1
	xn(1:nsd,1:nnc) = x_temp(1:nsd,offset:endset)

	if (myid.lt.numproc-1) then
	   call MPI_SEND(lock,1,MPI_DOUBLE_PRECISION,myid+1,101,
	1	MPI_COMM_WORLD,ierr)
	end if
	close(file)
	return
	end
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine readallx(xn)

        include "global.h"
        real* 8 xn(nsd,nn)
        integer lock,ierr,status(MPI_STATUS_SIZE)
        integer file,offset,endset,i
        real* 8 x_temp(nsd,nn),xn_temp(nsd,nnc)

        file=25
        open(file, FILE="mxyz.dat", STATUS="old")
        do i=1,nn
           read(file,*) xn(:,i)
        enddo
	close(file)
        return
        end

c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readien(ien)
	
	include "global.h"
	integer ien(nen,nec)
	integer lock,ierr,status(MPI_STATUS_SIZE)
	character*4 ifp
	integer file,offset,endset,ien_temp(nen,ne)

	file=21
	open(file, FILE="mien.dat", STATUS="old")
	do i=1,ne
	   read(file,*) ien_temp(:,i)
	enddo

	lock = 1
	if (myid.gt.0) then
	   call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1	MPI_COMM_WORLD,status,ierr)
	end if

	offset=myid*maxnec+1
	endset=offset+nec
	ien(:,1:nec)=ien_temp(:,offset:endset)

	if (myid.lt.numproc-1) then
	   call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1	MPI_COMM_WORLD,ierr)
	end if

	return
	end
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readrng(rngface)

	include "global.h"
	integer rngface(neface,nec), rngfaceieee(neface,nec/2+1)
	integer lock,ierr,io,status(MPI_STATUS_SIZE)
	character*4 ifp
	integer file,offset,endset,rng_temp(neface,ne)

	file=22
	open(file, FILE="mrng.dat", STATUS="old")
	
	do i=1,ne
	   read(file,*) rng_temp(:,i)
	enddo

c	maxrecl = neface * maxnec * 4
	lock = 1
	if (myid.gt.0) then
	   call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1	MPI_COMM_WORLD,status,ierr)
	end if

	offset=myid*maxnec+1
	endset=offset+nec
	rngface(:,1:nec)=rng_temp(:,offset:endset)

	if (myid.lt.numproc-1) then
	   call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1	MPI_COMM_WORLD,ierr)
	end if
	
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
	1	   MPI_COMM_WORLD,status,ierr)
	      mynrng = max(mynrng, nrng)
	   end do
	   nrng = mynrng
	   write(6,1000) nrng
 1000	   format("mrng: ",i2," boundaries")
	else
	   call MPI_SEND(mynrng,1,MPI_INTEGER,0,101,
	1	MPI_COMM_WORLD,ierr)
	   
	end if

        call MPI_BCAST(nrng,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	
	return
	end
