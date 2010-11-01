!c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Lucy Zhang
!c       Calculate deformation gradient for ALE.  7/1/99
!c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine defgrad(xloc, xrefloc,ien,f,jac,finv,jacinv)
	use fluid_variables, only: nen,ne,nn,nsd,nquad,wq,sq
	implicit none

	integer ien(nen,ne)
	real* 8 xloc(nsd,nn), xrefloc(nsd,nn)
	real* 8 x(nsd,nen), f(nsd,nsd,nquad,ne)
	real* 8 finv(nsd,nsd,nquad,ne)
	real* 8 eft0,det, jac(nquad,ne), minor(nsd), jacinv(nquad,ne)
	real* 8 sh(0:nsd,nen), ph(0:nsd,nen)
        real* 8 xrefx(nsd), xrefy(nsd), xrefz(nsd)
        real* 8 refxx(nsd), refxy(nsd), refxz(nsd)
	real* 8 xr(nsd,nsd), cf(nsd,nsd),sx(nsd,nsd)
	real* 8 f11,f12,f13,f21,f22,f23,f31,f32,f33
	real* 8 finv11,finv12,finv13,finv21,finv22,finv23,finv31
	real* 8 finv32,finv33
	integer inl, ie, isd, iq, node
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c.... calculate d(Na)/d(xref)
	do ie=1,ne
	   do inl=1,nen
	      do isd=1,nsd
		 x(isd,inl)=xrefloc(isd,ien(inl,ie))
	      enddo
	   enddo
	   do iq=1,nquad
	      do isd = 1,nsd
		 xrefx(isd) = 0.0
		 xrefy(isd) = 0.0
		 xrefz(isd) = 0.0
	      enddo
	      if (nen.eq.4) then
		 include "sh3d4n.h"
	      else if (nen.eq.8) then
		 include "sh3d8n.h"
	      end if

!c.... calculate the first derivative, d(x)/d(xref)
	      do inl=1,nen
		 do isd=1,nsd
		    xrefx(isd)=xrefx(isd)+sh(1,inl)*xloc(isd,ien(inl,ie))      
		    xrefy(isd)=xrefy(isd)+sh(2,inl)*xloc(isd,ien(inl,ie))      
		    xrefz(isd)=xrefz(isd)+sh(3,inl)*xloc(isd,ien(inl,ie))
		 enddo
	      enddo
	   
	      f(1,1,iq,ie) = xrefx(1)
	      f(1,2,iq,ie) = xrefy(1)
	      f(1,3,iq,ie) = xrefz(1)
	      f(2,1,iq,ie) = xrefx(2)
	      f(2,2,iq,ie) = xrefy(2)
	      f(2,3,iq,ie) = xrefz(2)
	      f(3,1,iq,ie) = xrefx(3)
	      f(3,2,iq,ie) = xrefy(3)
	      f(3,3,iq,ie) = xrefz(3)
	      f11 = f(1,1,iq,ie)
	      f12 = f(1,2,iq,ie)
	      f13 = f(1,3,iq,ie)
	      f21 = f(2,1,iq,ie)
	      f22 = f(2,2,iq,ie)
	      f23 = f(2,3,iq,ie)
	      f31 = f(3,1,iq,ie)
	      f32 = f(3,2,iq,ie)
	      f33 = f(3,3,iq,ie)
	      minor(1) = f22*f33 - f23*f32
	      minor(2) = f23*f31 - f21*f33
	      minor(3) = f21*f32 - f22*f31
	      jac(iq,ie) = f11*minor(1) + f12*minor(2) + f13*minor(3)
	      
	      finv(1,1,iq,ie) = (f22*f33-f23*f32)/jac(iq,ie)
	      finv(1,2,iq,ie) = (f23*f31-f21*f33)/jac(iq,ie)
	      finv(1,3,iq,ie) = (f21*f32-f22*f31)/jac(iq,ie)
	      finv(2,1,iq,ie) = (f13*f32-f12*f33)/jac(iq,ie)
	      finv(2,2,iq,ie) = (f11*f33-f13*f31)/jac(iq,ie)
	      finv(2,3,iq,ie) = (f12*f31-f11*f32)/jac(iq,ie)
	      finv(3,1,iq,ie) = (f12*f23-f13*f22)/jac(iq,ie)
	      finv(3,2,iq,ie) = (f13*f21-f11*f23)/jac(iq,ie)
	      finv(3,3,iq,ie) = (f11*f22-f12*f21)/jac(iq,ie)
	      finv11=finv(1,1,iq,ie)
	      finv12=finv(1,2,iq,ie)
	      finv13=finv(1,3,iq,ie)
	      finv21=finv(2,1,iq,ie)
	      finv22=finv(2,2,iq,ie)
	      finv23=finv(2,3,iq,ie)
	      finv31=finv(3,1,iq,ie)
	      finv32=finv(3,2,iq,ie)
	      finv33=finv(3,3,iq,ie)
	      minor(1) = finv22 * finv33 - finv23 * finv32 
	      minor(2) = finv23 * finv31 - finv21 * finv33 
	      minor(3) = finv21 * finv32 - finv22 * finv31 
	      jacinv(iq,ie)=finv11 *minor(1) &
		   +finv12 *minor(2)         &
		   +finv13 *minor(3)
	   enddo
	enddo
	return
	end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Calculate inverse of deformation gradient for ALE.  7/1/99
!c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine defgradinv(xloc,xrefloc,ien,finv,jacinv)
	use fluid_variables, only:nen,ne,nn,nsd,nquad,wq,sq
	implicit none

	integer ien(nen,ne)
	real* 8 xloc(nsd,nn), xrefloc(nsd,nn)
	real* 8 x(nsd,nen), finv(nsd,nsd,ne)

	real* 8 eft0,det, jac(ne), jacinv(ne),minor(nsd)
	real* 8 sh(0:nsd,nen), ph(0:nsd,nen)
        real* 8 refxx(nsd), refxy(nsd), refxz(nsd)
	real* 8 xr(nsd,nsd), cf(nsd,nsd),sx(nsd,nsd)
	integer inl, ie, isd, iq, node
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do ie=1,ne
!c.... calculate d(Na)/d(x)
           do inl=1,nen
              do isd=1,nsd
                 x(isd,inl)=xloc(isd,ien(inl,ie))
              enddo
           enddo
           do iq=1,nquad
	      do isd = 1,nsd
		 refxx(isd) = 0.0
		 refxy(isd) = 0.0
		 refxz(isd) = 0.0
	      enddo

	      if (nen.eq.4) then
		 include "sh3d4n.h"
	      else if (nen.eq.8) then
		 include "sh3d8n.h"
	      end if

!c.... calculate the first derivative, d(xref)/d(x)
	      do inl=1,nen
		 do isd=1,nsd
		    refxx(isd)=refxx(isd)+sh(1,inl)*xrefloc(isd,ien(inl,ie))  
		    refxy(isd)=refxy(isd)+sh(2,inl)*xrefloc(isd,ien(inl,ie))  
		    refxz(isd)=refxz(isd)+sh(3,inl)*xrefloc(isd,ien(inl,ie))
		 enddo
	      enddo
	      finv(1,1,ie) = refxx(1)
	      finv(1,2,ie) = refxy(1)
	      finv(1,3,ie) = refxz(1)
	      finv(2,1,ie) = refxx(2)
	      finv(2,2,ie) = refxy(2)
	      finv(2,3,ie) = refxz(2)
	      finv(3,1,ie) = refxx(3)
	      finv(3,2,ie) = refxy(3)
	      finv(3,3,ie) = refxz(3)
	      minor(1) = finv(2,2,ie)*finv(3,3,ie)-finv(2,3,ie)*finv(3,2,ie)
	      minor(2) = finv(2,3,ie)*finv(3,1,ie)-finv(2,1,ie)*finv(3,3,ie)
	      minor(3) = finv(2,1,ie)*finv(3,2,ie)-finv(2,2,ie)*finv(3,1,ie)
	      jacinv(ie)=finv(1,1,ie)*minor(1) &
		   +finv(1,2,ie)*minor(2)      &
		   +finv(1,3,ie)*minor(3)
	   enddo
	enddo

       return
       end
