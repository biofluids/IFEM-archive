!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Lucy Zhang
!   Calculate mesh velocities for ALE.  7/1/99
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine velocity(xn,xnold,meshvel,refvel,finv,dloc,ien)
  use fluid_variables
  use run_variables, only: dt
  use global_constants
  implicit none

  real(8),intent(in) :: xn(nsd,nn),xnold(nsd,nn)
  real(8),intent(out):: meshvel(nsd,nn),refvel(nsd,nquad,ne)
  real(8),intent(in) :: finv(nsd,nsd,nquad,ne)
  real(8),intent(in) :: dloc(ndf,nn)
  integer,intent(in):: ien(nen,ne)
  
  !real(8) meshveln(nsd,nn),
  real(8) :: dx
  real(8) :: convel(nsd)
  !real(8) :: f(nsd,nsd,nquad,ne)
  
  real(8) :: d(ndf,nen),mv(nsd,nen),diffv(nsd,nen)
  integer :: inn, isd,i,j,ie,inl,iq
  real(8) :: eft0,sh(0:nsdpad,nenpad)
  real(8) :: xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)
  real(8) :: x(nsd,nen),xloc(nsd,nn),det
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!... calculate mesh velocities, vhat.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  refvel(1:nsd,1:nquad,1:ne) = 0.0d0
  do inn = 1,nn
     do isd = 1,nsd
        dx = xn(isd,inn) - xnold(isd,inn)
        meshvel(isd,inn) = dx/dt
     enddo
     
  enddo
  !write(*,'("meshvel ",3F12.8)') meshvel(1:3,66)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!... calculate referential velocities, w
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
        if (nen.eq.4) then
           include "sh3d4n.h"
        elseif (nen.eq.8) then
           include "sh3d8n.h"
        endif
        eft0=abs(det)*wq(iq)

        do isd=1,nsd
           convel(isd)=0.0
        enddo
        do inl=1,nen
           convel(xsd)=convel(xsd)+sh(0,inl)*diffv(xsd,inl)
           convel(ysd)=convel(ysd)+sh(0,inl)*diffv(ysd,inl)
           convel(zsd)=convel(zsd)+sh(0,inl)*diffv(zsd,inl)
        enddo

        refvel(1:nsd,iq,ie)=0
        
        do i = 1,nsd
           refvel(i,iq,ie) = 0.0d0
           do j = 1,nsd
              refvel(i,iq,ie) = refvel(i,iq,ie)+finv(i,j,iq,ie)*convel(j)
           enddo
        enddo

     enddo
  enddo

  return
end subroutine velocity





