!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  S. Aliabadi                                                          c
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine vol(xloc, ien) 
  use fluid_variables
  implicit none


  integer ien(nen,ne)
  real* 8 xloc(nsd,nn)
  real* 8 x(nsdpad,nenpad)

  real* 8 eft0,det
  real* 8 sh(0:nsdpad,nenpad)
  real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

!  real* 8 e_gas,e_liq,p_gas,p_liq
  real* 8 e_liq, p_liq
  integer inl,ie,iq,isd

  !integer ir,status(MPI_STATUS_SIZE)

  p_liq =  0.0

  do ie=1,ne

     do inl=1,nen
        do isd=1,nsd
           x(isd,inl) = xloc(isd,ien(inl,ie))
	    enddo
	 enddo

	 e_liq = 0.0
	 do iq=1,nquad
	    if (nen.eq.4) then
	       include "sh3d4n.h"
	    else if (nen.eq.8) then
		   include "sh3d8n.h"
	    end if
	    eft0 = abs(det) * wq(iq)  

	    e_liq = e_liq + eft0

	 enddo

	 p_liq = p_liq + e_liq
	   
  enddo
  
  liq=p_liq
	
  return
end subroutine vol
