c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine formd(xn,ds,rngface,ien,hn,hm)

      implicit none
	  include "global.h"
      include "malloc.h"

	  integer  rngface(neface,nec),ien(nen,nec)
      real* 8  xn(nsd,nnc),ds(ndf,nnc)
	  real* 8  hn(nnc),hm(nn_loc)
	  integer idf, inl, iec, irng, ieface, inface, inn
      logical assemble
      real* 8 eps1,eps2

      real* 8 hs(nrng,nnc), h(nrng,nn_loc)
      pointer (hsptr, hs),(hptr, h)

c  construct lien and mien mappings from ien array
	  hsptr = malloc(nrng*nnc*fsize)
	  hptr  = malloc(nrng*nn_loc*fsize)

	  eps1 = -1000000.0 
	  eps2 = -10000.0 

	  call fclear (h,nrng*nn_loc)

	  do inn=1,nnc
		do idf=1,ndf
		  ds(idf,inn) = eps1
		enddo
	  enddo

	  do ieface=1,neface
		do inface=1,nnface
		  inl = map(ieface,inface,etype)
		  do iec=1,nec
			irng = rngface(ieface,iec)
			if(irng.ne.0) then
			  h(irng,ien(inl,iec)) = h(irng,ien(inl,iec)) + 1.0                
			endif
		  enddo
		enddo
	  enddo

	  assemble=.true.
	  call scatter(h,hs,nrng,assemble,hn,hm)

	  do irng=1,nrng
        do  inn=1,nnc
		  do  idf=1,ndf
			if((hs(irng,inn).gt.1.0e-8).and.(bc(idf,irng).gt.0)) 
	1			 ds(idf,inn) = bv(idf,irng)
		  enddo
        enddo
	  enddo

	  do inn=1,nnc
        do idf=1,ndf
		  if(ds(idf,inn).lt.eps2) ds(idf,inn) = ic(idf)
        enddo
	  enddo

	  if(static) then
        do inn = 1,nnc
		  do idf = 1,nsd 
			ds(idf,inn) = 0.0
		  enddo 
        enddo
	  endif

	  call free(hsptr)
	  call free(hptr)

	  return
	  end
c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine formid(ids, rngface, ien, hn,hm)

      implicit none
	  include "global.h"
      include "malloc.h"

	  integer ids(ndf,nnc),rngface(neface,nec),ien(nen,nec)
	  real* 8  hn(nnc),hm(nn_loc)
	  integer idf, inl, iec, irng, ieface, inface, inn, iflag
      logical assemble
      real* 8 epsr,epsl
      integer ierr,status(MPI_STATUS_SIZE)

      real* 8 ds(ndf,nnc),d(ndf,nn_loc)
      pointer (dsptr, ds),(dptr, d)

      dsptr = malloc(ndf*nnc*fsize)
      dptr  = malloc(ndf*nn_loc*fsize)

      call fclear (d,ndf*nn_loc)
      epsr = 0.0001         
      epsl = 0.000001      

	  do ieface=1,neface
		do inface=1,nnface
		  inl = map(ieface,inface,etype)
		  do iec=1,nec
			irng = rngface(ieface,iec)
			if(irng.ne.0) then
			  do idf = 1,ndf
				if(d(idf,ien(inl,iec)).lt.epsr) 
	1				 d(idf,ien(inl,iec)) = bc(idf,irng)+epsl
			  enddo
			endif
		  enddo
		end do
	  end do

	  assemble=.true.
	  call scatter(d,ds,ndf,assemble,hn,hm)

	  do inn=1,nnc
		do idf=1,ndf
		  ids(idf,inn) = ds(idf,inn)
		enddo
	  enddo

	  call free(dsptr)
	  call free(dptr)

	  if(static) then
        do inn = 1,nnc
		  do idf = 1,nsd 
			ids(idf,inn) = 1
		  enddo 
        enddo
	  endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  nq = 0
	  iflag = 0
	  if(myid.eq.0) iflag = 1

 888  continue
	  if(iflag.eq.1) then
		do inn=1,nnc
		  do idf=1,ndf
			if(ids(idf,inn).eq.0) then
			  nq = nq + 1
			  ids(idf,inn) = nq
			else
			  ids(idf,inn) = 0
			endif
		  enddo
        enddo
        if(myid.ne.numproc-1) call MPI_SEND(nq,1,MPI_INTEGER,myid+1,
	1		 101,MPI_COMM_WORLD,ierr)
        goto 999  
	  else 
        call MPI_RECV(nq,1,MPI_INTEGER,myid-1,101,
	1		 MPI_COMM_WORLD,status,ierr)
        iflag = 1
        goto 888
	  endif
 999  continue

	  call MPI_BCAST(nq,1,MPI_INTEGER,numproc-1,MPI_COMM_WORLD,ierr)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	  return
	  end