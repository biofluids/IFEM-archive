      subroutine writelu

c  For processors with id < ndf, write lu decomposition data
c  for future use: amat, icoln, ikeep
      include "global.h"

      integer idf, iidf
      integer lockrec,ierr,status(MPI_STATUS_SIZE)

      if (myid.lt.ndf) then
        idf = myid+1
                    
c  Write amat
        lockrec = 0
        if (myid.gt.0) then
          call MPI_RECV(lockrec,1,MPI_INTEGER,myid-1,101,
     &         MPI_COMM_WORLD,status,ierr)
        end if
        if (ecounttot(idf).gt.0) then
          call ewd_open("mamat", ifp)
          call ewd_lseek(ifp, lockrec*8, 0)
          call ewd_write(ifp, amat, licn*8)
          call ewd_close(ifp)
          lockrec = lockrec + licn
        end if
        if (myid.lt.ndf-1) then
          call MPI_SEND(lockrec,1,MPI_INTEGER,myid+1,101,
     &         MPI_COMM_WORLD,ierr)
        end if
        
c  Write icoln
        lockrec = 0
        if (myid.gt.0) then
          call MPI_RECV(lockrec,1,MPI_INTEGER,myid-1,101,
     &         MPI_COMM_WORLD,status,ierr)
        end if
        if (ecounttot(idf).gt.0) then
          call ewd_open("micol", ifp)
          call ewd_lseek(ifp, lockrec*8, 0)
          call ewd_write(ifp, icoln, licn*8)
          call ewd_close(ifp)
          lockrec = lockrec + licn
        end if
        if (myid.lt.ndf-1) then
          call MPI_SEND(lockrec,1,MPI_INTEGER,myid+1,101,
     &         MPI_COMM_WORLD,ierr)
        end if
        
c  Write ikeep
        lockrec = 0
        if (myid.gt.0) then
          call MPI_RECV(lockrec,1,MPI_INTEGER,myid-1,101,
     &         MPI_COMM_WORLD,status,ierr)
        end if
        if (ecounttot(idf).gt.0) then
          call ewd_open("mikeep", ifp)
          call ewd_lseek(ifp, lockrec*8, 0)
          call ewd_write(ifp, ikeep, 5*order*8)
          call ewd_close(ifp)
          lockrec = lockrec + 5*order
        end if
        if (myid.lt.ndf-1) then
          call MPI_SEND(lockrec,1,MPI_INTEGER,myid+1,101,
     &         MPI_COMM_WORLD,ierr)
        end if

      end if

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine readlu

c  For processors with id < ndf, read lu decomposition data
c  from files: amat, icoln, ikeep
      include "global.h"

      integer idf, iidf
      integer lockrec,ierr,status(MPI_STATUS_SIZE)

      if (myid.lt.ndf) then
        idf = myid+1
                    
c  Read amat
        lockrec = 0
        if (myid.gt.0) then
          call MPI_RECV(lockrec,1,MPI_INTEGER,myid-1,101,
     &         MPI_COMM_WORLD,status,ierr)
        end if
        if (ecounttot(idf).gt.0) then
          call ewd_open("mamat", ifp)
          call ewd_lseek(ifp, lockrec*8, 0)
          call ewd_read(ifp, amat, licn*8)
          call ewd_close(ifp)
          lockrec = lockrec + licn
        end if
        if (myid.lt.ndf-1) then
          call MPI_SEND(lockrec,1,MPI_INTEGER,myid+1,101,
     &         MPI_COMM_WORLD,ierr)
        end if
        
c  Read icoln
        lockrec = 0
        if (myid.gt.0) then
          call MPI_RECV(lockrec,1,MPI_INTEGER,myid-1,101,
     &         MPI_COMM_WORLD,status,ierr)
        end if
        if (ecounttot(idf).gt.0) then
          call ewd_open("micol", ifp)
          call ewd_lseek(ifp, lockrec*8, 0)
          call ewd_read(ifp, icoln, licn*8)
          call ewd_close(ifp)
          lockrec = lockrec + licn
        end if
        if (myid.lt.ndf-1) then
          call MPI_SEND(lockrec,1,MPI_INTEGER,myid+1,101,
     &         MPI_COMM_WORLD,ierr)
        end if
        
c  Read ikeep
        lockrec = 0
        if (myid.gt.0) then
          call MPI_RECV(lockrec,1,MPI_INTEGER,myid-1,101,
     &         MPI_COMM_WORLD,status,ierr)
        end if
        if (ecounttot(idf).gt.0) then
          call ewd_open("mikeep", ifp)
          call ewd_lseek(ifp, lockrec*8, 0)
          call ewd_read(ifp, ikeep, 5*order*8)
          call ewd_close(ifp)
          lockrec = lockrec + 5*order
        end if
        if (myid.lt.ndf-1) then
          call MPI_SEND(lockrec,1,MPI_INTEGER,myid+1,101,
     &         MPI_COMM_WORLD,ierr)
        end if

      end if

      return
      end


        
