	subroutine diskin(x,d,f,dd)
	include "global.h"
	real* 8 x(nsd,nnc),d(ndf,nnc),f(nnc),dd(ndf+4,nnc)
      integer ndft,i,j
      integer lk,ir,status(MPI_STATUS_SIZE)

ccccccccccccc
c..... ndf+1= # of deg. of freedom+phi; ndf+3= # of dof+coordinates
      ndft = ndf + 4
	call fclear (dd,ndft*nnc)
	maxrecl = ndft* maxnnc * 8
	lk = 1
	if (myid.gt.0) then
	 call MPI_RECV(lk,1,MPI_INTEGER,myid-1,5,MPI_COMM_WORLD,status,ir)
	end if
	call ewd_open("data",ifp)
	call ewd_lseek(ifp, myid*maxrecl, 0)
	call ewd_read(ifp,dd, ndft*nnc*8)
	call ewd_close(ifp)
	if (myid.lt.numproc-1) then
	   call MPI_SEND(lk,1,MPI_INTEGER,myid+1,5,MPI_COMM_WORLD,ir)
	end if
ccccccccccccc
	do i=1,nnc
	   do j=1,nsd
	      x(j,i)=dd(j,i)
	   enddo
	   do j=1,ndf
	      d(j,i)=dd(j+3,i)
	   enddo
	   f(i)=dd(ndft,i)
	enddo
	f(i) = dd(ndft,i)

ccccccccccccc

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine diskout(x,d,a,dd)

	include "global.h"
	real* 8 x(nsd,nnc),d(ndf,nnc),a(nnc),dd(ndf+4,nnc)
	character*9 dout
	integer i,j,ndft,i1,i2,i3,i4
	integer lk,ir,status(MPI_STATUS_SIZE)

	ndft = ndf + 4
ccccccccccccc
	do i=1,nnc
	   do j=1,nsd
	      dd(j,i)=x(j,i)
	      
	   enddo
	   do j=1,ndf
	      dd(j+3,i) = d(j,i) 
	   enddo	
	   dd(ndft,i) = a(i) 
	enddo
ccccccccccccc

      ibase = ichar('0')    !! integer value for char '0'

	dout = "data.0000"

      i4 = idisk/1000
      i3 = (idisk-i4*1000)/100
      i2 = (idisk-i4*1000-i3*100)/10
      i1 = (idisk-i4*1000-i3*100-i2*10)/1
      
      i4 = i4 + ibase
      i3 = i3 + ibase
      i2 = i2 + ibase
      i1 = i1 + ibase

      dout(6:6) = char(i4)
      dout(7:7) = char(i3)
      dout(8:8) = char(i2)
      dout(9:9) = char(i1)
ccccccccc
	maxrecl = ndft* maxnnc * 8
      lk = 1
      if (myid.gt.0) then
      call MPI_RECV(lk,1,MPI_INTEGER,myid-1,5,MPI_COMM_WORLD,status,ir)
	end if
	call ewd_open(dout, ifp)
	call ewd_lseek(ifp, myid*maxrecl, 0)
	call ewd_write(ifp,dd,ndft*nnc*8)
	call ewd_close(ifp)
      if (myid.lt.numproc-1) then
      call MPI_SEND(lk,1,MPI_INTEGER,myid+1,5,MPI_COMM_WORLD,ir)
	end if
ccccccccc

	return
	end



	subroutine diskout2(d)

	include "global.h"
	real* 8 d(nsd,nnc)
	character*9 dout
	integer i,j,i1,i2,i3,i4
	integer lk,ir,status(MPI_STATUS_SIZE)

cccccccccccccccccccccccccc

	ibase = ichar('0')	!! integer value for char '0'

	dout = "mesh.0000"

	i4 = idisk/1000
	i3 = (idisk-i4*1000)/100
	i2 = (idisk-i4*1000-i3*100)/10
	i1 = (idisk-i4*1000-i3*100-i2*10)/1
      
	i4 = i4 + ibase
	i3 = i3 + ibase
	i2 = i2 + ibase
	i1 = i1 + ibase

	dout(6:6) = char(i4)
	dout(7:7) = char(i3)
	dout(8:8) = char(i2)
	dout(9:9) = char(i1)
ccccccccc
	maxrecl = nsd* maxnnc * 8
	lk = 1
	if (myid.gt.0) then
	   call MPI_RECV(lk,1,MPI_INTEGER,myid-1,5,MPI_COMM_WORLD,status,ir)
	end if
	
	call ewd_open(dout, ifp)
	write(*,*) 'in mesh.f ifp=',dout
	call ewd_lseek(ifp, myid*maxrecl, 0)
	call ewd_write(ifp,d,nsd*nnc*8)
	call ewd_close(ifp)

	if (myid.lt.numproc-1) then
	   call MPI_SEND(lk,1,MPI_INTEGER,myid+1,5,MPI_COMM_WORLD,ir)
	end if
ccccccccc

	return
	end
