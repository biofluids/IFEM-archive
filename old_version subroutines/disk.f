	subroutine diskin(d,dd)
	include "global.h"
	real* 8 d(ndf,nn),dd(ndf,nn)
	integer ndft,i,j
	integer lk,ir
c	integer status(MPI_STATUS_SIZE)
	
ccccccccccccc
	ndft = ndf
	call fclear (dd,ndft*nn)
	maxrecl = ndft* maxnn * 8
	lk = 1

	call readrestart(myid*maxrecl,0,dd, ndft*nn*8)
	call readrestart(0,0,dd, ndft*nn*8)

ccccccccccccc
	do i=1,nn
	   do j=1,ndf
	      d(j,i) = dd(j,i)
	   enddo
	enddo
ccccccccccccc
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine diskout(d,dd)
	
	include "global.h"
	real* 8 d(ndf,nn),dd(ndf,nn)
	character*9 dout
	integer i,j,ndft,i1,i2,i3,i4
	integer lk,ir
c	integer status(MPI_STATUS_SIZE)
	
	ndft = ndf
ccccccccccccc
	do i=1,nn
	   do j=1,ndf
	      dd(j,i) = d(j,i) 
	   enddo
	enddo
ccccccccccccc
	
	ibase = ichar('0')	!! integer value for char '0'
	
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
	maxrecl = ndft* maxnn * 8
	lk = 1

c	call datawrite(i4,i3,i2,i1,myid*maxrecl,0,dd,ndft*nn*8)
	call datawrite(i4,i3,i2,i1,0,0,dd,ndft*nn*8)

ccccccccc
	return
	end
