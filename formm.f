c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine formdm(xn,ds,rngface,ien,diameter)

	implicit none
	include "global.h"
	include "malloc.h"

	integer  rngface(neface,ne),ien(nen,ne)
	real* 8  xn(nsd,nn),ds(nsd,nn),dsf(nn)
	real* 8  tt
	integer isd, inl, iec, irng, ieface, inface, inn
	logical assemble
	real* 8 eps1,eps2,diameter,x0,w0,damp,wd
	real* 8 h(nrng,nn)
	pointer (hptr, h)

c  construct lien and mien mappings from ien array

        hptr  = malloc(nrng*nn*fsize)

        eps1 = -1000000.0 
        eps2 = -10000.0 
	w0=sqrt(spring/mass)
	damp=damper/mass/2/w0
	wd=w0*sqrt(1-damp**2)
	x0=amp*diameter
	call fclear (h,nrng*nn)

        do inn=1,nn
	   do isd=1,nsd
	      ds(isd,inn) = eps1
	   enddo
	enddo

	do ieface=1,neface
	   do inface=1,nnface
	      inl = map(ieface,inface,etype)
	      do iec=1,ne
		 irng = rngface(ieface,iec)
		 if(irng.ne.0) then
		    h(irng,ien(inl,iec)) = h(irng,ien(inl,iec)) + 1.0  
		 endif
	      enddo
	   enddo
	enddo

	do irng=1,nrng
	  do inn=1,nn
	     do  isd=1,nsd
c..... the boundary value, bvd is inputted for oscillation frequency.
c.....  the displacement is a periodic sine wave, a function of time.
	      if((h(irng,inn).gt.1.0e-8).and.(bcd(isd,irng).gt.0)) then
		     ds(isd,inn)=x0*sin(bvd(isd,irng)*tt)
c		     ds(isd,inn)=exp(-w0*damp*tt)*(x0*cos(wd*tt)+w0*damp
c	1	      *x0/wd*sin(wd*tt))
		   if (bvd(isd,irng).eq.0.0)  ds(isd,inn)=0.0
	      endif

	     enddo
	  enddo
	enddo
	do inn = 1,nn
	   do isd = 1,nsd 
		 if(ds(isd,inn).lt.eps1) ds(isd,inn) = 0.0
	   enddo 
	enddo

        call free(hptr)

	return
	end
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine formidm(ids, rngface, ien)

	implicit none
	include "global.h"
	include "malloc.h"

	integer ids(nsd,nn),rngface(neface,ne),ien(nen,ne)
	integer isd, inl, iec, irng, ieface, inface, inn, iflag
	logical assemble
	real* 8 epsr,epsl
	integer ierr,status(MPI_STATUS_SIZE)

	real* 8 d(nsd,nn)
	pointer (dptr, d)

	real* 8 dsf(nn)
	pointer (dsfptr, dsf)

	dptr  = malloc(nsd*nn*fsize)
	dsfptr= malloc(nn*fsize)

	call fclear (d,nsd*nn)
	epsr = 0.0001         
	epsl = 0.000001      
	do ieface=1,neface
	   do iec=1,ne
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

	do inn=1,nn
	   do isd=1,nsd
	      ids(isd,inn) = d(isd,inn)
	   enddo
	enddo

        call free(dptr)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      nq = 0

	do inn=1,nn
	   do isd=1,nsd
		 if(ids(isd,inn).eq.0) then
		    nq = nq + 1
		    ids(isd,inn) = nq
		 else
		    ids(isd,inn) = 0
		 endif
	   enddo
	enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	return
	end