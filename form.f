c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine formd(xn,ds,rngface,ien)

	implicit none
	include "global.h"
	include "malloc.h"

	integer  rngface(neface,ne),ien(nen,ne)
	real* 8  xn(nsd,nn),ds(ndf,nn)
	integer idf, inl, iec, irng, ieface, inface, inn
c	logical assemble
	real* 8 eps1,eps2

	real* 8 hs(nrng,nn), h(nrng,nn)
	pointer (hsptr, hs),(hptr, h)

c  construct lien and mien mappings from ien array
	hsptr = malloc(nrng*nn*fsize)
	hptr  = malloc(nrng*nn*fsize)

	eps1 = -1000000.0 
	eps2 = -10000.0 

	call fclear (h,nrng*nn)

	do inn=1,nn
	   do idf=1,ndf
	      ds(idf,inn) = eps1
	   enddo
	enddo

	do ieface=1,neface
	   do inface=1,nnface
	      inl = map(ieface,inface,etype)
	      do iec=1,ne
			irng = rngface(ieface,iec)
			if (irng.ne.0) then
				h(irng,ien(inl,iec)) = h(irng,ien(inl,iec)) + 1.0
			endif
	      enddo
	   enddo
	enddo

	call equal(h,hs,nrng*nn)

	do irng=1,nrng
	   do inn=1,nn
	      do idf=1,ndf
		  if((hs(irng,inn).gt.1.0e-8).and.(bc(idf,irng).gt.0)) 
	1	      ds(idf,inn) = bv(idf,irng)
	      enddo
	   enddo
	enddo

	do inn=1,nn
	   do idf=1,ndf
	      if(ds(idf,inn).lt.eps2) ds(idf,inn) = ic(idf)
	   enddo
	enddo

	if(static) then
	   do inn = 1,nn
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
	subroutine formid(ids, rngface, ien)

	implicit none
	include "global.h"
	include "malloc.h"

	integer ids(ndf,nn),rngface(neface,ne),ien(nen,ne)
	integer idf, inl, iec, irng, ieface, inface, inn, iflag
	logical assemble
	real* 8 epsr,epsl

	real* 8 ds(ndf,nn),d(ndf,nn)
	pointer (dsptr, ds),(dptr, d)

	dsptr = malloc(ndf*nn*fsize)
	dptr  = malloc(ndf*nn*fsize)

	call fclear (d,ndf*nn)
	epsr = 0.0001         
	epsl = 0.000001      

	do ieface=1,neface
	   do inface=1,nnface
	      inl = map(ieface,inface,etype)
	      do iec=1,ne
		    irng = rngface(ieface,iec)
		    if(irng.ne.0) then
		      do idf = 1,ndf
		        if(d(idf,ien(inl,iec)).lt.epsr) 
	1		      d(idf,ien(inl,iec)) = bc(idf,irng)+epsl
		      enddo
		    endif
	      enddo
	   enddo
	enddo

	call equal(d,ds,ndf*nn)

	do inn=1,nn
	   do idf=1,ndf
	      ids(idf,inn) = ds(idf,inn)
	   enddo
	enddo

	call free(dsptr)
	call free(dptr)
	
	if(static) then
	   do inn = 1,nn
	      do idf = 1,nsd 
		    ids(idf,inn) = 1
	      enddo 
	   enddo
	endif
	
	nq = 0

	do inn=1,nn
	   do idf=1,ndf
	      if(ids(idf,inn).eq.0) then
		    nq = nq + 1
		    ids(idf,inn) = nq
	      else
			ids(idf,inn) = 0
	      endif
	   enddo
	enddo
	
	return
	end
