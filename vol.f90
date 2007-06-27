!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  S. Aliabadi                                                          c
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine vol(xloc, ien) 
  use fluid_variables
  implicit none

  integer :: ien(nen,ne)
  real(8) :: xloc(nsd,nn)
  real(8) :: x(nsd,nen)
  real(8) :: eft0,det
  real(8) :: sh(0:nsd,nen)
  real(8) :: xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)
  real(8) :: e_liq,p_liq
  integer :: inl,ie,iq

  p_liq =  0.0

  do ie=1,ne
     do inl=1,nen
        x(:,inl) = xloc(:,ien(inl,ie))
     enddo

     e_liq = 0.0
     do iq=1,nquad
	 	if (nsd==2) then
		    if (nen.eq.3) then !calculate shape function at quad point
			   include "sh2d3n.h"
			elseif (nen.eq.4) then
				include "sh2d4n.h"
			endif
		elseif (nsd==3) then
		    if (nen.eq.4) then !calculate shape function at quad point
			   include "sh3d4n.h"
			elseif (nen.eq.8) then
				include "sh3d8n.h"
			endif
		endif

        eft0 = abs(det) * wq(iq)  

        e_liq = e_liq + eft0

     enddo

     p_liq = p_liq + e_liq
         
  enddo
  
  liq=p_liq
  
       
  return
end subroutine vol
