c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine formdm(xn,ds,rngface,ien,hn,hm,diameter)

	implicit none
	include "global.h"
	include "malloc.h"

	integer  rngface(neface,nec),ien(nen,nec)
	real* 8  xn(nsd,nnc),ds(nsd,nnc),dsf(nnc)
	real* 8  hn(nnc),hm(nn_loc),tt
	integer isd, inl, iec, irng, ieface, inface, inn
	logical assemble
	real* 8 eps1,eps2,diameter,x0,w0,damp,wd
	real* 8 hs(nrng,nnc), h(nrng,nn_loc)
	pointer (hsptr, hs),(hptr, h)

c  construct lien and mien mappings from ien array
        hsptr = malloc(nrng*nnc*fsize)
        hptr  = malloc(nrng*nn_loc*fsize)

        eps1 = -1000000.0 
        eps2 = -10000.0 
	w0=sqrt(spring/mass)
	damp=damper/mass/2/w0
	wd=w0*sqrt(1-damp**2)
	x0=amp*diameter

	call fclear (h,nrng*nn_loc)

c	write(*,*) 'amp=',amp
        do inn=1,nnc
	   do isd=1,nsd
	      ds(isd,inn) = eps1
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
	  do inn=1,nnc
	     do  isd=1,nsd
c..... the boundary value, bvd is inputted for oscillation frequency.
c.....  the displacement is a periodic sine wave, a function of time.
	      if((hs(irng,inn).gt.1.0e-8).and.(bcd(isd,irng).gt.0)) then
		     ds(isd,inn)=x0*sin(bvd(isd,irng)*tt)
c		     ds(isd,inn)=exp(-w0*damp*tt)*(x0*cos(wd*tt)+w0*damp
c	1	      *x0/wd*sin(wd*tt))
		   if (bvd(isd,irng).eq.0.0)  ds(isd,inn)=0.0
	      endif

	     enddo
	  enddo
	enddo
	do inn = 1,nnc
	   do isd = 1,nsd 
		 if(ds(isd,inn).lt.eps1) ds(isd,inn) = 0.0
	   enddo 
	enddo
        call free(hsptr)
        call free(hptr)

	return
	end
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine formidm(ids, rngface, ien, hn,hm)

	implicit none
	include "global.h"
	include "malloc.h"

	integer ids(nsd,nnc),rngface(neface,nec),ien(nen,nec)
	real* 8  hn(nnc),hm(nn_loc)
	integer isd, inl, iec, irng, ieface, inface, inn, iflag
	logical assemble
	real* 8 epsr,epsl
	integer ierr,status(MPI_STATUS_SIZE)

	real* 8 ds(nsd,nnc),d(nsd,nn_loc)
	pointer (dsptr, ds),(dptr, d)

	real* 8 dsf(nnc),df(nn_loc)
	pointer (dsfptr, dsf),(dfptr, df)

	dsptr = malloc(nsd*nnc*fsize)
	dptr  = malloc(nsd*nn_loc*fsize)
	dsfptr= malloc(nnc*fsize)
	dfptr = malloc(nn_loc*fsize)

	call fclear (d,nsd*nn_loc)
	epsr = 0.0001         
	epsl = 0.000001      
	do ieface=1,neface
	   do iec=1,nec
	      irng = rngface(ieface,iec)
	      if(irng.ne.0) then
		 do inface=1,nnface
		    inl = map(ieface,inface,etype)
		    do isd = 1,nsd
		       if(d(isd,ien(inl,iec)).lt.epsr) 
	1		    d(isd,ien(inl,iec)) = bcd(isd,irng)+epsl
		    enddo
		 enddo
	      endif
	   enddo
	enddo
        assemble=.true.
        call scatter(d,ds,nsd,assemble,hn,hm)

	do inn=1,nnc
	   do isd=1,nsd
	      ids(isd,inn) = ds(isd,inn)
	   enddo
	enddo
        call free(dsptr)
        call free(dptr)
        call free(dsfptr)
        call free(dfptr)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        nq = 0
        iflag = 0
        if(myid.eq.0) iflag = 1

 888    continue
        if(iflag.eq.1) then
	   do inn=1,nnc
	      do isd=1,nsd
		 if(ids(isd,inn).eq.0) then
		    nq = nq + 1
		    ids(isd,inn) = nq
		 else
		    ids(isd,inn) = 0
		 endif
	      enddo
	   enddo
	   if(myid.ne.numproc-1) call MPI_SEND(nq,1,MPI_INTEGER,myid+1,
	1	101,MPI_COMM_WORLD,ierr)
	   goto 999  
        else 
	   call MPI_RECV(nq,1,MPI_INTEGER,myid-1,101,
	1	MPI_COMM_WORLD,status,ierr)
	   iflag = 1
	   goto 888
        endif
 999    continue

        call MPI_BCAST(nq,1,MPI_INTEGER,numproc-1,MPI_COMM_WORLD,ierr)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	return
	end



