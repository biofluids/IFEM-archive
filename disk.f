	  subroutine diskin(d,dd)
	  include "global.h"
	  real* 8 d(ndf,nnc),dd(ndf,nnc)
      integer ndft,i,j
      integer lk,ir,status(MPI_STATUS_SIZE)

ccccccccccccc
      ndft = ndf
	  call fclear (dd,ndft*nnc)
	  maxrecl = ndft* maxnnc * 8
      lk = 1
      if (myid.gt.0) then
		call MPI_RECV(lk,1,MPI_INTEGER,myid-1,5,MPI_COMM_WORLD,status,ir)
	  end if
	  call ewd_open("dddd",ifp)
	  call ewd_lseek(ifp, myid*maxrecl, 0)
	  call ewd_read(ifp,dd, ndft*nnc*8)
	  call ewd_close(ifp)
      if (myid.lt.numproc-1) then
		call MPI_SEND(lk,1,MPI_INTEGER,myid+1,5,MPI_COMM_WORLD,ir)
	  end if
ccccccccccccc
      do i=1,nnc
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
	  real* 8 d(ndf,nnc),dd(ndf,nnc)
      character*9 dout
      integer i,j,ndft,i1,i2,i3,i4
      integer lk,ir,status(MPI_STATUS_SIZE)

      ndft = ndf
ccccccccccccc
      do i=1,nnc
		do j=1,ndf
		  dd(j,i) = d(j,i) 
		enddo
	  enddo
ccccccccccccc

      ibase = ichar('0')		!! integer value for char '0'

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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine restartout(d,dd)

	  include "global.h"
	  real* 8 d(ndf,nnc),dd(ndf,nnc)
      character*9 dout
      integer i,j,ndft,i1,i2,i3,i4
      integer lk,ir,status(MPI_STATUS_SIZE)

      ndft = ndf
ccccccccccccc
      do i=1,nnc
		do j=1,ndf
		  dd(j,i) = d(j,i) 
		enddo
	  enddo
ccccccccccccc

      ibase = ichar('0')		!! integer value for char '0'

	  dout = "dddd"

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
