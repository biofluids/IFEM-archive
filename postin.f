      subroutine postin(x,d,f,dd,file)
      include "global.h"
      real* 8 x(nsd,nnc),dd(ndf+4,nnc)
      real* 8 d(ndf,nnc),f(nnc)
      character file*(*)
      integer ndft,i,j
      integer lk,ir,status(MPI_STATUS_SIZE)

      ndft=ndf+4
      call fclear (dd,ndft*nnc)

      maxrecl = ndft* maxnnc * 8
      lk = 1
      if (myid.gt.0) then
        call MPI_RECV(lk,1,MPI_INTEGER,myid-1,5,MPI_COMM_WORLD,status,ir)
      end if

      call ewd_open(file,ifp)
      call ewd_lseek(ifp, myid*maxrecl, 0)
      call ewd_read(ifp,dd, ndft*nnc*8)
      call ewd_close(ifp)
      if (myid.lt.numproc-1) then
        call MPI_SEND(lk,1,MPI_INTEGER,myid+1,5,MPI_COMM_WORLD,ir)
      end if
      do i=1,nnc
         do j=1,nsd
c            x(j,i)=dd(j,i)
            dd(j,i)=x(j,i)+dd(j,i)
         enddo
         f(i)=dd(ndft,i)
         do j=1,ndf
            d(j,i)=dd(j+3,i)
         enddo
      enddo
	return
	end


      subroutine postin2(d,file)
      include "global.h"
      real* 8 d(nsd,nnc),f(nnc)
      character file*(*)
      integer ndf,i,j
      integer lk,ir,status(MPI_STATUS_SIZE)

c	write(*,*) 'file=',file
      call fclear (d,nsd*nnc)
      maxrecl = nsd* maxnnc * 8
      lk = 1
      if (myid.gt.0) then
        call MPI_RECV(lk,1,MPI_INTEGER,myid-1,5,MPI_COMM_WORLD,status,ir)
      end if

      call ewd_open(file,ifp)
      call ewd_lseek(ifp, myid*maxrecl, 0)
      call ewd_read(ifp,d, nsd*nnc*8)
      call ewd_close(ifp)
      if (myid.lt.numproc-1) then
        call MPI_SEND(lk,1,MPI_INTEGER,myid+1,5,MPI_COMM_WORLD,ir)
      end if

	return
	end




