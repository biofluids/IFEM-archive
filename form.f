c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine formd(xn,ds,dsf,rngface,ien)

      implicit none
	include "global.h"
      include "malloc.h"

	integer  rngface(neface,ne),ien(nen,ne)
      real* 8  xn(nsd,nn),ds(ndf,nn),dsf(nn)
	integer idf, inl, iec, irng, ieface, inface, inn
      real* 8 eps1,eps2

      real* 8 h(nrng,nn)
      pointer (hptr, h)

c       construct lien and mien mappings from ien array

        hsptr = malloc(nrng*nn*fsize)
        hptr  = malloc(nrng*nn*fsize)

        eps1 = -1000000.0 
        eps2 = -10000.0 

	  call fclear (h,nrng*nn)

      do inn=1,nn
	  do idf=1,ndf
		ds(idf,inn) = eps1
	  enddo
        dsf (inn) = eps1
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
        do  inn=1,nn
          do  idf=1,ndf
		  if((h(irng,inn).gt.1.0e-8).and.(bc(idf,irng).gt.0)) 
     &                             ds(idf,inn) = bv(idf,irng)
          enddo
          if((h(irng,inn).gt.1.0e-8).and.(bcf(irng).gt.0))
     &                             dsf(inn) = bvf(irng)
        enddo
      enddo

      do inn=1,nn
        do idf=1,ndf
          if(ds(idf,inn).lt.eps2) ds(idf,inn) = ic(idf)
        enddo
        if(dsf(inn).lt.eps2) dsf(inn) = icf
      enddo

        do inn=1,nn
           if(xn(1,inn).lt.interface(1)) dsf(inn) = 1.0
           if(xn(2,inn).lt.interface(2)) dsf(inn) = 1.0
           if(xn(3,inn).lt.interface(3)) dsf(inn) = 1.0
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
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine formid(ids, idsf, rngface, ien)

      implicit none
	include "global.h"
      include "malloc.h"

	integer ids(ndf,nn),idsf(nn),rngface(neface,ne),ien(nen,ne)
	integer idf, inl, iec, irng, ieface, inface, inn, iflag
      logical assemble
      real* 8 epsr,epsl

      real* 8 d(ndf,nn)
      pointer (dptr, d)

      real* 8 df(nn)
      pointer (dfptr, df)

      dptr  = malloc(ndf*nn*fsize)
      dfptr = malloc(nn*fsize)

      call fclear (d,ndf*nn)
	call fclear (df,nn)
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
     &		     d(idf,ien(inl,iec)) = bc(idf,irng)+epsl
			enddo
	  	    if(df(ien(inl,iec)).lt.epsr) df(ien(inl,iec)) = bcf(irng)+epsl
		  endif
		enddo
	  enddo
	enddo

      do inn=1,nn
	  do idf=1,ndf
	    ids(idf,inn) = d(idf,inn)
	  enddo
        idsf(inn)= df(inn)
	enddo

        call free(dptr)
        call free(dfptr)

        if(static) then
        do inn = 1,nn
        do idf = 1,nsd 
        ids(idf,inn) = 1
        enddo 
        enddo
        endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        nqf= 0

	do inn=1,nn
	  if(idsf(inn).eq.0) then
		nqf = nqf + 1
		idsf(inn) = nqf
	  else
	  idsf(inn) = 0
	  endif
      enddo
        
	return
	end
