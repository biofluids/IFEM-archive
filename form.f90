module form
  implicit none


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formd(xn,ds,rngface,ien)
  use fluid_variables

  implicit none

  !include "malloc.fi"

  integer  rngface(neface,ne),ien(nen,ne)
  real* 8  xn(nsd,nn),ds(ndf,nn)
  integer idf, inl, iec, irng, ieface, inface, inn
!	logical assemble
  real* 8 eps1,eps2

  real*8 :: hs(nrng,nn), h(nrng,nn)
!  pointer (hsptr, hs),(hptr, h)

 !...construct lien and mien mappings from ien array
!  hsptr = malloc(nrng*nn*fsize)
!  hptr  = malloc(nrng*nn*fsize)

  eps1 = -1000000.0 
  eps2 = -10000.0 

  h(:,:) = 0.0d0

  do inn=1,nn
     do idf=1,ndf
        ds(idf,inn) = eps1
     enddo
  enddo

  do ieface=1,neface
     do inface=1,nnface
        inl = mapping(ieface,inface,etype)
        do iec=1,ne
           irng = rngface(ieface,iec)
		   if (irng.ne.0) then
		      h(irng,ien(inl,iec)) = h(irng,ien(inl,iec)) + 1.0
		   endif
	    enddo
	 enddo
  enddo

  !call equal(h,hs,nrng*nn)
  hs(:,:) = h(:,:)

  do irng=1,nrng
     do inn=1,nn
        do idf=1,ndf
           if (hs(irng,inn).gt.1.0e-8) then
              if (bc(idf,irng) .gt. 0) then
                 ds(idf,inn) = bv(idf,irng)
              endif
	       endif
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

!  call free(hsptr)
!  call free(hptr)

  return
end subroutine formd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formid(ids, rngface, ien)
  use fluid_variables
  implicit none

  !include "malloc.fi"

  integer ids(ndf,nn),rngface(neface,ne),ien(nen,ne)
  integer idf, inl, iec, irng, ieface, inface, inn
!  logical assemble
  real* 8 epsr,epsl

  real* 8 ds(ndf,nn),d(ndf,nn)
!  pointer (dsptr, ds),(dptr, d)

!  dsptr = malloc(ndf*nn*fsize)
!  dptr  = malloc(ndf*nn*fsize)

  d(1:ndf,1:nn) = 0.0d0

  epsr = 0.0001         
  epsl = 0.000001      

  do ieface=1,neface
     do inface=1,nnface
        inl = mapping(ieface,inface,etype)
        do iec=1,ne
           irng = rngface(ieface,iec)
		   if(irng.ne.0) then
		      do idf = 1,ndf
		         if(d(idf,ien(inl,iec)).lt.epsr) d(idf,ien(inl,iec)) = bc(idf,irng)+epsl
		      enddo
		   endif
	    enddo
	 enddo
  enddo

  !call equal(d,ds,ndf*nn)
  ds(:,:) = d(:,:)

  do inn=1,nn
     do idf=1,ndf
        ids(idf,inn) = ds(idf,inn)
     enddo
  enddo

!  call free(dsptr)
!  call free(dptr)
	
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
end subroutine formid



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formbc(ds,rngface,ien)
  use global_constants
  use run_variables, only: tt
  use fluid_variables
  implicit none
  !include "malloc.fi"

	integer  rngface(neface,ne),ien(nen,ne)
	real*8  xn(nsd,nn),ds(ndf,nn)
	integer idf, inl, iec, irng, ieface, inface, inn
!	logical assemble
	real*8 eps1,eps2,tt_ramp

	real*8 hs(nrng,nn), h(nrng,nn)

	!real*8,parameter :: pi=3.14159265d0

 !...construct lien and mien mappings from ien array


  

  eps1 = -1000000.0 
  eps2 = -10000.0 

  h(1:nrng,1:nn) = 0.0d0
  
	!do inn=1,nn
	!   do idf=1,ndf
	!      ds(idf,inn) = eps1
	!   enddo
	!enddo

	do ieface=1,neface
	   do inface=1,nnface
	      inl = mapping(ieface,inface,etype)
	      do iec=1,ne
			irng = rngface(ieface,iec)
			if (irng.ne.0) then
				h(irng,ien(inl,iec)) = h(irng,ien(inl,iec)) + 1.0
			endif
	      enddo
	   enddo
	enddo

	!call equal(h,hs,nrng*nn)
	hs(:,:) = h(:,:)


    tt_ramp = 0.05
    if ((tt >= 0).and.(tt <= tt_ramp)) then
         bv(udf,2) = 0.5  * tt/tt_ramp
		 bv(udf,4) =-0.5  * tt/tt_ramp  !...inlet (6) velocity in direction #degree of freedom#
	elseif ((tt >= tt_ramp).and.(tt <= (2*tt_ramp))) then
         bv(udf,2) = 0.5
		 bv(udf,4) =-0.5
    endif

	!write(*,'("boundary 6 (inflow) x-velocity: ",f7.4," cm/s")')bv(udf,6)
	!write(*,'("boundary 6 (inflow) z-velocity: ",f6.3," cm/s")')bv(wdf,6)


	do irng=1,nrng
	   do inn=1,nn
	      do idf=1,ndf
		     if (hs(irng,inn).gt.1.0e-8) then
                  if (bc(idf,irng) .gt. 0) then
                     ds(idf,inn) = bv(idf,irng)
                  endif
	         endif
	      enddo
	   enddo
	enddo

	!do inn=1,nn
	!   do idf=1,ndf
	!      if(ds(idf,inn).lt.eps2) ds(idf,inn) = ic(idf)
	!write(*,*) "attention (formbc): low values for d"
	!   enddo
	!enddo

	if(static) then
       ds(1:ndf,1:nn) = 0.0
	endif

  return
end subroutine formbc

end module form


