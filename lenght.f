c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	S. Aliabadi                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine lenght(xloc,ien,hg) 

      implicit none
	include "global.h"

      integer ien(nen,nec)
	real* 8 xloc(nsd,nn_loc),x(nsdpad,nenpad)
	real* 8 hg(nec)

      real* 8 eft0,det
      real* 8 sh(0:nsdpad,nenpad)
      real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

	real* 8 evol
	real* 8 gmin,gmax,lmin,lmax
	integer i,inl,ie,iq,isd

      integer ierr,status(MPI_STATUS_SIZE)

      vmin = +10000000.0
	vmax = -10000000.0
      hmin = +10000000.0
	hmax = -10000000.0

      do ie=1,nec

        do inl=1,nen
        do isd=1,nsd
        x(isd,inl) = xloc(isd,ien(inl,ie))
        enddo
        enddo

        evol = 0.0
        do iq=1,nquad
	  if (nen.eq.4) then
		include "sh3d4n.h"
	  else if (nen.eq.8) then
		include "sh3d8n.h"
	  end if

        eft0 = abs(det) * wq(iq)  
	  evol = evol + eft0
        enddo

	  vmin = min(vmin,evol)
	  vmax = max(vmax,evol)

	  if(hg_vol) then
	  hg(ie) = evol**(1.0/3.0)
	  if(twod) hg(ie) = sqrt(evol)
	  if(nen.eq.4) hg(ie) = (8.0*evol)**(1.0/3.0)
	  else
	  call get_hg(x, hg(ie))
	  endif

	  hg(ie) = delta(0)*hg(ie)

	  hmin = min(hmin,hg(ie))
	  hmax = max(hmax,hg(ie))

        enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      lmin = hmin
      lmax = hmax

	gmin = lmin
	gmax = lmax
	if (myid.eq.0) then
	do i=1,numproc-1
      call MPI_RECV(lmax,1,MPI_DOUBLE_PRECISION,i,101,
     &                     MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(lmin,1,MPI_DOUBLE_PRECISION,i,102,
     &                     MPI_COMM_WORLD,status,ierr)
	gmax = max(gmax,lmax)
	gmin = min(gmin,lmin)
	end do
	else
      call MPI_SEND(lmax,1,MPI_DOUBLE_PRECISION,0,101,
     &                     MPI_COMM_WORLD,ierr)
      call MPI_SEND(lmin,1,MPI_DOUBLE_PRECISION,0,102,
     &                     MPI_COMM_WORLD,ierr)
	endif


      call MPI_BCAST(gmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(gmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	hmin = gmin
	hmax = gmax
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      lmin = vmin
      lmax = vmax

	gmin = lmin
	gmax = lmax
	if (myid.eq.0) then
	do i=1,numproc-1
      call MPI_RECV(lmax,1,MPI_DOUBLE_PRECISION,i,101,
     &                     MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(lmin,1,MPI_DOUBLE_PRECISION,i,102,
     &                     MPI_COMM_WORLD,status,ierr)
	gmax = max(gmax,lmax)
	gmin = min(gmin,lmin)
	end do
	else
      call MPI_SEND(lmax,1,MPI_DOUBLE_PRECISION,0,101,
     &                     MPI_COMM_WORLD,ierr)
      call MPI_SEND(lmin,1,MPI_DOUBLE_PRECISION,0,102,
     &                     MPI_COMM_WORLD,ierr)
	endif


      call MPI_BCAST(gmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(gmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	vmin = gmin
	vmax = gmax
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      return
      end
