      subroutine postin(d,f,dd,file)
      include "global.h"
      real* 8 d(ndf,nnc),f(nnc),dd(ndf+1,nnc)
      character file*(*)
      integer ndft,i,j
      integer lk,ir,status(MPI_STATUS_SIZE)

      ndft = ndf + 1
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
        f(i) = dd(ndft,i)
        do j=1,ndf
          d(j,i) = dd(j,i)
        enddo
      enddo
      
	return
	end
