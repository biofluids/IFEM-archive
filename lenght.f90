!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  S. Aliabadi                                                          c
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine lenght(xloc,ien,hg) 
  use fluid_variables
  implicit none

  integer :: ien(nen,ne)
  real(8) :: xloc(nsd,nn),x(nsdpad,nenpad)
  real(8) :: hg(ne)

  real(8) :: eft0,det
  real(8) :: sh(0:nsd,nen)
  real(8) :: xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

  real(8) :: evol
  real(8) :: gmin,gmax,lmin,lmax
  integer :: inl,ie,iq,isd

  !integer ierr,status(MPI_STATUS_SIZE)

  vmin = +10000000.0
  vmax = -10000000.0
  hmin = +10000000.0
  hmax = -10000000.0

  do ie=1,ne
  
     do inl=1,nen
        do isd=1,nsd
           x(isd,inl) = xloc(isd,ien(inl,ie))
        enddo
     enddo

     evol = 0.0
     do iq=1,nquad
	   if (nsd==2) then
		 if (nen==3) then 
		    include "sh2d3n.h"
		 elseif (nen==4) then
		    include "sh2d4n.h"
	     endif
	   elseif (nsd==3) then
         if (nen.eq.4) then
           include "sh3d4n.h"
         else if (nen.eq.8) then
           include "sh3d8n.h"
         end if
	   endif

        eft0 = abs(det) * wq(iq)  
        evol = evol + eft0
     enddo

     vmin = min(vmin,evol)
     vmax = max(vmax,evol)

     if(hg_vol) then
        if (nsd==3) then
			hg(ie) = evol**(1.0/3.0)
			if(nen == 4) hg(ie) = (8.0*evol)**(1.0/3.0)
		elseif (nsd==2) then
		    hg(ie) = sqrt(evol)
		endif
     else
        call get_hg(x, hg(ie))
     endif

     hg(ie) = delta(0)*hg(ie)

     hmin = min(hmin,hg(ie))
     hmax = max(hmax,hg(ie))

  enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  lmin = hmin
  lmax = hmax

  gmin = lmin
  gmax = lmax
!c  if (myid.eq.0) then
!c     do i=1,numproc-1
!c        call MPI_RECV(lmax,1,MPI_DOUBLE_PRECISION,i,101,
!c  1      MPI_COMM_WORLD,status,ierr)
!c        call MPI_RECV(lmin,1,MPI_DOUBLE_PRECISION,i,102,
!c  1      MPI_COMM_WORLD,status,ierr)
          gmax = max(gmax,lmax)
          gmin = min(gmin,lmin)
!c	   end do
!c	else
!c	   call MPI_SEND(lmax,1,MPI_DOUBLE_PRECISION,0,101,
!c	1	MPI_COMM_WORLD,ierr)
!c	   call MPI_SEND(lmin,1,MPI_DOUBLE_PRECISION,0,102,
!c	1	MPI_COMM_WORLD,ierr)
!c	endif


!c	call MPI_BCAST(gmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!c	call MPI_BCAST(gmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!c	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  hmin = gmin
  hmax = gmax
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  lmin = vmin
  lmax = vmax

  gmin = lmin
  gmax = lmax
!c	if (myid.eq.0) then
!c	   do i=1,numproc-1
!c	      call MPI_RECV(lmax,1,MPI_DOUBLE_PRECISION,i,101,
!c	1	   MPI_COMM_WORLD,status,ierr)
!c	      call MPI_RECV(lmin,1,MPI_DOUBLE_PRECISION,i,102,
!c	1	   MPI_COMM_WORLD,status,ierr)
          gmax = max(gmax,lmax)
          gmin = min(gmin,lmin)
!c	   end do
!c	else
!c	   call MPI_SEND(lmax,1,MPI_DOUBLE_PRECISION,0,101,
!c	1	MPI_COMM_WORLD,ierr)
!c	   call MPI_SEND(lmin,1,MPI_DOUBLE_PRECISION,0,102,
!c	1	MPI_COMM_WORLD,ierr)
!c	endif


!c	call MPI_BCAST(gmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!c	call MPI_BCAST(gmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!c	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  vmin = gmin
  vmax = gmax
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  return
end subroutine lenght
