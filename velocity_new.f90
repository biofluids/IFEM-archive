!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Lucy Zhang
!c   Calculate mesh velocities for ALE.  7/1/99
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine velocity_new(xn,xnold,meshvel,refvel, &
          dloc,ien)
      use fluid_variables, only: nsd,nen,ne,nn,nquad,wq,sq,ndf
      use global_constants
      use run_variables, only: dt
      implicit none
      
      integer ien(nen,ne)
      real* 8 xn(nsd,nn), xnold(nsd,nn)
      real* 8 meshveln(nsd,nn),dx,refvel(nsd,nquad,ne)
      real* 8 meshvel(nsd,nn),convel(nsd)
      !real* 8 f(nsd,nsd,nquad,ne),finv(nsd,nsd,nquad,ne)
      real* 8 d(ndf,nen),mv(nsd,nen),dloc(ndf,nn),diffv(nsd,nen)
      integer inn, isd,node,i,j,ie,inl,iq
      real* 8 eft0,sh(0:nsd,nen)
      real* 8 xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)
      real* 8 x(nsd,nen),xloc(nsd,nn),det
      real(8) tmp_xw
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c... calculate mesh velocities, vhat.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call fclear(refvel,nsd*nquad*ne)
      do inn = 1,nn
         do isd = 1,nsd
            dx = xn(isd,inn) - xnold(isd,inn)
            meshvel(isd,inn) = dx/dt
         enddo
!c        write(*,*) 'meshveln=',meshveln(1,inn),meshveln(2,inn),meshveln(3,inn)
      enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c... calculate referential velocities, w
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do ie = 1,ne
         do inl=1,nen
            do isd=1,nsd
               d(isd,inl)=dloc(isd,ien(inl,ie))
               x(isd,inl)=xloc(isd,ien(inl,ie))
               mv(isd,inl)=meshvel(isd,ien(inl,ie))
               diffv(isd,inl)=d(isd,inl)-mv(isd,inl)
            enddo
         enddo
         do iq=1,nquad
	
	if (nsd == 3) then  ! 3-D case
            if (nen.eq.4) then
               include "sh3d4n.h"
            elseif (nen.eq.8) then
               include "sh3d8n.h"
            endif
	end if

	if (nsd == 2) then ! 2-D case
	   if (nen.eq.3) then
                  include "sh2d3n.h"
              elseif (nen.eq.4) then
                  include "sh2d4n.h"
              endif   
	end if

            eft0=abs(det)*wq(iq)

            do isd=1,nsd
               convel(isd)=0.0
            enddo

            do inl=1,nen
		do isd=1,nsd
	           convel(isd)=convel(isd)+sh(0,inl)*diffv(isd,inl)
		end do
!               convel(xsd)=convel(xsd)+sh(0,inl)*diffv(xsd,inl)
!               convel(ysd)=convel(ysd)+sh(0,inl)*diffv(ysd,inl)
!               convel(zsd)=convel(zsd)+sh(0,inl)*diffv(zsd,inl)
            enddo

            refvel(isd,iq,ie)=0
        
            do i = 1,nsd
               do j = 1,nsd
		!tmp_xw=finv(i,j,iq,ie)*convel(j)
		tmp_xw=convel(j)
                refvel(i,iq,ie)=refvel(i,iq,ie)+tmp_xw
               enddo
            enddo
         enddo
      enddo

      return
      end





