!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readx(xn)
      use fluid_variables, only: nsd,nn

	real* 8 xn(nsd,nn)
!	integer lock,ierr,status(MPI_STATUS_SIZE)
	integer file,offset,endset,i
!	real* 8 x_temp(nsd,nn),xn_temp(nsd,nn)

	file=23
	open(file, FILE="mxyz.dat", STATUS="old")

	do i=1,nn
	   read(file,*) xn(1:3,i)
	enddo	
	close(file)
	return
	end
!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readien(ien)
	use fluid_variables, only: nen,ne

	integer ien(nen,ne)
!	integer lock,ierr,status(MPI_STATUS_SIZE)
!	character*4 ifp
	integer file
	!integer offset,endset
	!integer :: ien_temp(nen,ne)

	file=21
	open(file, FILE="mien.dat", STATUS="old")
	do i=1,ne
	   read(file,*) ien(:,i)
	enddo
	close(file)
	return
	end
!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readrng(rngface)
	use fluid_variables, only: ne,neface

	integer rngface(neface,ne)
	!integer rngfaceieee(neface,ne/2+1)
!	integer lock,ierr,io,status(MPI_STATUS_SIZE)
	!character*4 ifp
	integer file
	!integer offset,endset

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
