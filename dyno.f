	subroutine dyno(xloc,dloc,rng,ien)

	implicit none
	include "global.h"

      integer seven
      parameter(seven = 7)
	real*8 tmp(seven),tot(seven)

	integer rng(neface,nec),ien(nen,nec)       
	real*8 xloc(nsd,nn_loc),dloc(ndf,nn_loc)

	real*8 x(nsdpad,10),d(10)
	real*8 xr(nsdpad,nsdpad) 
	real*8 pp
	real*8 n1,n2,n3,c1,c2,c3
	real*8 dfx,dfy,dfz,dmx,dmy,dmz,dea
	integer i,inl,ie,isd,ieface,inface,isurf

	integer ierr,io,status(MPI_STATUS_SIZE)

      call fclear(tmp,seven) 

      do ie=1,nec
	do ieface = 1,neface
	do isurf = 1,surf(0)
	if(rng(ieface,ie).eq.surf(isurf)) then
	  do inface=1,nnface
	  inl = map(ieface,inface,etype)
        do isd=1,nsd
        x(isd,inface) = xloc(isd,ien(inl,ie))
        enddo
        d(inface) = dloc(pdf,ien(inl,ie))
        enddo

	  c1 = 0.0
	  c2 = 0.0
	  c3 = 0.0
	  pp = 0.0
	  do inface=1,nnface
	  c1 = c1 + x(1,inface)
	  c2 = c2 + x(2,inface)
	  c3 = c3 + x(3,inface)
	  pp = pp + d(  inface)
	  enddo
	  c1 = c1/nnface
	  c2 = c2/nnface
	  c3 = c3/nnface
	  pp = pp/nnface

	  if (nen.eq.4) then
	  xr(1,1)=x(1,1)-x(1,3)
	  xr(2,1)=x(2,1)-x(2,3)
	  xr(3,1)=x(3,1)-x(3,3)
	  xr(1,2)=x(1,2)-x(1,3)
	  xr(2,2)=x(2,2)-x(2,3)
	  xr(3,2)=x(3,2)-x(3,3)
	  n1=0.5*(-xr(2,1)*xr(3,2)+xr(2,2)*xr(3,1))
	  n2=0.5*(-xr(3,1)*xr(1,2)+xr(1,1)*xr(3,2))
	  n3=0.5*(-xr(1,1)*xr(2,2)+xr(1,2)*xr(2,1))
	  else
	  xr(1,1)=0.5*(x(1,2)+x(1,3)-x(1,1)-x(1,4))
	  xr(2,1)=0.5*(x(2,2)+x(2,3)-x(2,1)-x(2,4))
	  xr(3,1)=0.5*(x(3,2)+x(3,3)-x(3,1)-x(3,4))
	  xr(1,2)=0.5*(x(1,3)+x(1,4)-x(1,1)-x(1,2))
	  xr(2,2)=0.5*(x(2,3)+x(2,4)-x(2,1)-x(2,2))
	  xr(3,2)=0.5*(x(3,3)+x(3,4)-x(3,1)-x(3,2))
	  n1=(-xr(2,1)*xr(3,2)+xr(2,2)*xr(3,1))
	  n2=(-xr(3,1)*xr(1,2)+xr(1,1)*xr(3,2))
	  n3=(-xr(1,1)*xr(2,2)+xr(1,2)*xr(2,1))
	  endif

	  dfx = -pp*n1
	  dfy = -pp*n2
	  dfz = -pp*n3
	  dmx = dfz*c2-dfy*c3
	  dmy = dfx*c3-dfz*c1
	  dmz = dfy*c1-dfx*c2
	  dea = sqrt(n1**2+n2**2+n3**2)

        tmp(1) = tmp(1)+dfx
        tmp(2) = tmp(2)+dfy
        tmp(3) = tmp(3)+dfz
        tmp(4) = tmp(4)+dmx
	  tmp(5) = tmp(5)+dmy
	  tmp(6) = tmp(6)+dmz
	  tmp(7) = tmp(7)+dea
      endif
      enddo
      enddo
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
	do i = 1,seven
	tot(i) = tmp(i)
	enddo
      if (myid.eq.0) then
        do io=1,numproc-1
             call MPI_RECV(tmp,seven,MPI_DOUBLE_PRECISION,io,
     &                         101,MPI_COMM_WORLD,status,ierr)
	       do i = 1,seven
	         tot(i) = tot(i)+tmp(i)
	       enddo
        enddo
        do io=1,numproc-1
             call MPI_SEND(tot,seven,MPI_DOUBLE_PRECISION,io,
     &                         102,MPI_COMM_WORLD,ierr)
        enddo
        else
             call MPI_SEND(tmp,seven,MPI_DOUBLE_PRECISION,0,
     &                         101,MPI_COMM_WORLD,ierr)
             call MPI_RECV(tot,seven,MPI_DOUBLE_PRECISION,0,
     &                         102,MPI_COMM_WORLD,status,ierr)
      endif

        if(myid.eq.0) then
        write(9,9) tt,tot(7),tot(1),tot(2),tot(3)
        write(8,9) tt,tot(7),tot(4),tot(5),tot(6)
        endif
   9    format(5(2x,e12.5))

	return
	end
