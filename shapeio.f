c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine readtmp(shrk,shrknode,cnn,ncnn,cnn2,ncnn2)
	  
	  include "global.h"
	  real* 8 shrk(0:nsd,maxconn,nquad*nec)
	  real* 8 shrknode(maxconn,nnc)
	  integer cnn(maxconn,nqdc), ncnn(nqdc) 
	  integer cnn2(maxconn,nnc), ncnn2(nnc) 
	  
      integer lock,ierr,status(MPI_STATUS_SIZE)

	  maxrecl = (nsd+1)*maxconn * nquad*maxnec * 8
      lock = 1
	  if (myid.gt.0) then
		call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
	  end if
	  call ewd_open("mshp", ifp)
	  call ewd_lseek(ifp, myid*maxrecl,0)
	  call ewd_read(ifp, shrk,(nsd+1)*maxconn*nquad*nec*8)
	  call ewd_close(ifp)
	  if (myid.lt.numproc-1) then
		call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1		 MPI_COMM_WORLD,ierr)
	  end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  maxrecl = maxconn * maxnnc * 8
      lock = 1
	  if (myid.gt.0) then
		call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
	  end if
	  call ewd_open("mshn", ifp)
	  call ewd_lseek(ifp,myid*maxrecl,0)
	  call ewd_read(ifp,shrknode,maxconn*nnc*8)
	  call ewd_close(ifp)
	  if (myid.lt.numproc-1) then
		call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1		 MPI_COMM_WORLD,ierr)
	  end if
cccccccccccccccccccccccccccccccccccccccccccccccc
	  maxrecl = maxconn*maxnqdc*8
      lock = 1
	  if (myid.gt.0) then
		call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
	  end if
	  call ewd_open("mcnn", ifp)
	  call ewd_lseek(ifp, myid*maxrecl,0)
	  call ewd_read(ifp, cnn,maxconn*nqdc*8)
	  call ewd_close(ifp)
	  if (myid.lt.numproc-1) then
		call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1		 MPI_COMM_WORLD,ierr)
	  end if

cccccccccccccccccccccccccccccccccccccccccccccccc
	  maxrecl = maxnqdc*8
      lock = 1
	  if (myid.gt.0) then
		call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
	  end if
	  call ewd_open("mncc", ifp)
	  call ewd_lseek(ifp, myid*maxrecl,0)
	  call ewd_read(ifp,ncnn,nqdc*8)
	  call ewd_close(ifp)
	  if (myid.lt.numproc-1) then
		call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1		 MPI_COMM_WORLD,ierr)
	  end if
cccccccccccccccccccccccccccccccccccccccccccccccc
	  maxrecl = maxconn*maxnnc*8
      lock = 1
	  if (myid.gt.0) then
		call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
	  end if
	  call ewd_open("mcnn2", ifp)
	  call ewd_lseek(ifp, myid*maxrecl,0)
	  call ewd_read(ifp, cnn2,maxconn*nnc*8)
	  call ewd_close(ifp)
	  if (myid.lt.numproc-1) then
		call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1		 MPI_COMM_WORLD,ierr)
	  end if

cccccccccccccccccccccccccccccccccccccccccccccccc
	  maxrecl = maxnnc*8
      lock = 1
	  if (myid.gt.0) then
		call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
	  end if
	  call ewd_open("mncc2", ifp)
	  call ewd_lseek(ifp, myid*maxrecl,0)
	  call ewd_read(ifp,ncnn2,nnc*8)
	  call ewd_close(ifp)
	  if (myid.lt.numproc-1) then
		call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1		 MPI_COMM_WORLD,ierr)
	  end if

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

	  subroutine writetmp(shrk,shrknode,cnn,ncnn,cnn2,ncnn2)

	  include "global.h"
	  real* 8 shrk(0:nsd,maxconn,nquad*nec)
	  real* 8 shrknode(maxconn,nnc)
	  integer cnn(maxconn,nqdc), ncnn(nqdc) 
	  integer cnn2(maxconn,nnc), ncnn2(nnc) 
	  
      integer lock,ierr,status(MPI_STATUS_SIZE)

	  maxrecl = (nsd+1)*maxconn * nquad*maxnec * 8
      lock = 1
	  if (myid.gt.0) then
		call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
	  end if
	  call ewd_open("mshp", ifp)
	  call ewd_lseek(ifp, myid*maxrecl,0)
	  call ewd_write(ifp, shrk,(nsd+1)*maxconn*nquad*nec*8)
	  call ewd_close(ifp)
	  if (myid.lt.numproc-1) then
		call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1		 MPI_COMM_WORLD,ierr)
	  end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  maxrecl = maxconn * maxnnc * 8
      lock = 1
	  if (myid.gt.0) then
		call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
	  end if
	  call ewd_open("mshn", ifp)
	  call ewd_lseek(ifp, myid*maxrecl,0)
	  call ewd_write(ifp, shrknode,maxconn*nnc*8)
	  call ewd_close(ifp)
	  if (myid.lt.numproc-1) then
		call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1		 MPI_COMM_WORLD,ierr)
	  end if
cccccccccccccccccccccccccccccccccccccccccccccccc
	  maxrecl = maxconn*maxnqdc*8
      lock = 1
	  if (myid.gt.0) then
		call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
	  end if
	  call ewd_open("mcnn", ifp)
	  call ewd_lseek(ifp, myid*maxrecl,0)
	  call ewd_write(ifp, cnn,maxconn*nqdc*8)
	  call ewd_close(ifp)
	  if (myid.lt.numproc-1) then
		call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1		 MPI_COMM_WORLD,ierr)
	  end if

cccccccccccccccccccccccccccccccccccccccccccccccc
	  maxrecl = maxnqdc*8
      lock = 1
	  if (myid.gt.0) then
		call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
	  end if
	  call ewd_open("mncc", ifp)
	  call ewd_lseek(ifp, myid*maxrecl,0)
	  call ewd_write(ifp,ncnn,nqdc*8)
	  call ewd_close(ifp)
	  if (myid.lt.numproc-1) then
		call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1		 MPI_COMM_WORLD,ierr)
	  end if
cccccccccccccccccccccccccccccccccccccccccccccccc
	  maxrecl = maxconn*maxnnc*8
      lock = 1
	  if (myid.gt.0) then
		call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
	  end if
	  call ewd_open("mcnn2", ifp)
	  call ewd_lseek(ifp, myid*maxrecl,0)
	  call ewd_write(ifp, cnn2,maxconn*nnc*8)
	  call ewd_close(ifp)
	  if (myid.lt.numproc-1) then
		call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1		 MPI_COMM_WORLD,ierr)
	  end if

cccccccccccccccccccccccccccccccccccccccccccccccc
	  maxrecl = maxnnc*8
      lock = 1
	  if (myid.gt.0) then
		call MPI_RECV(lock,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
	  end if
	  call ewd_open("mncc2", ifp)
	  call ewd_lseek(ifp, myid*maxrecl,0)
	  call ewd_write(ifp,ncnn2,nnc*8)
	  call ewd_close(ifp)
	  if (myid.lt.numproc-1) then
		call MPI_SEND(lock,1,MPI_INTEGER,myid+1,101,
	1		 MPI_COMM_WORLD,ierr)
	  end if

      return
      end
