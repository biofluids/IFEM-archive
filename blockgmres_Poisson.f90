!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!residual for Laplace eqn!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  blockgmres_Poisson(cal_ne,cal_domain,xloc,I_fluid_den,p_inter,ien_den)
  use global_constants
  use run_variables
  use fluid_variables,only:nsd,iquad
!  use centermesh_variables
  use denmesh_variables, only:nn_den,ne_den,nen_den
  implicit none

  integer ien_den(nen_den,ne_den)
  real(8) xloc(nsd,nn_den)
  real(8) I_fluid_den(nn_den)
  real(8) p_inter(nn_den)

  real(8) x(nsd,nen_den),d(nen_den)
  real(8) eft0,det,effd,effm,effc
  real(8) sh(0:nsd,nen_den),ph(0:nsd,nen_den)
  real(8) xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real(8) dr(nsd)
  real(8) hg
  integer cal_ne
  integer cal_domain(cal_ne)
  
  integer inl,ie,isd,node,iq,ine

  integer,parameter :: ndfpad=5,nsdpad=3,nenpad=8,nquadpad=8
  integer nquad
  real(8) sq(0:nsdpad,nenpad,nquadpad)
  real(8) xq(nsdpad,nquadpad),wq(nquadpad)

  if(nsd==3) then
     write(*,*)'no 3d in block laplace'
     stop
  end if

!  if (nen_center==3) then
!     call quad2d3n(iquad,nquad,xq,wq,nsdpad,nquadpad)
!  end if
!  do iq=1,nquad
!     if(nen_center==3) then
!        sq(0,1,iq) = xq(1,iq)
!        sq(0,2,iq) = xq(2,iq)
!        sq(0,3,iq) = 1 - xq(1,iq) - xq(2,iq)
!     end if
!  end do
  if (nen_den==3) then
       call quad2d3n(iquad, nquad, xq, wq, nsdpad, nquadpad)
  else if (nen_den==4) then
       call quad2d4n(iquad, nquad, xq, wq, nsdpad, nquadpad)
  end if
  do iq=1,nquad
       if(nen_den==3) then
                  sq(0,1,iq) = xq(1,iq)
                  sq(0,2,iq) = xq(2,iq)
                  sq(0,3,iq) = 1 - xq(1,iq) - xq(2,iq)
        elseif (nen_den==4) then
                  sq(0,1,iq) = (1 - xq(1,iq)) * (1 - xq(2,iq)) / 4
                  sq(0,2,iq) = (1 + xq(1,iq)) * (1 - xq(2,iq)) / 4
                  sq(0,3,iq) = (1 + xq(1,iq)) * (1 + xq(2,iq)) / 4
                  sq(0,4,iq) = (1 - xq(1,iq)) * (1 + xq(2,iq)) / 4
                
                  sq(1,1,iq) = - (1 - xq(2,iq)) / 4
                  sq(1,2,iq) = + (1 - xq(2,iq)) / 4
                  sq(1,3,iq) = + (1 + xq(2,iq)) / 4
                  sq(1,4,iq) = - (1 + xq(2,iq)) / 4
        
                  sq(2,1,iq) = - (1 - xq(1,iq)) / 4
                  sq(2,2,iq) = - (1 + xq(1,iq)) / 4
                  sq(2,3,iq) = + (1 + xq(1,iq)) / 4
                  sq(2,4,iq) = + (1 - xq(1,iq)) / 4

        endif
  enddo

!==============================================================
!  do ie=1,ne_den     !loop over elements
  do ine=1,cal_ne
     ie=cal_domain(ine)

     do inl=1,nen_den
	x(1:nsd,inl) = xloc(1:nsd,ien_den(inl,ie))
	d(inl) = I_fluid_den(ien_den(inl,ie))
     end do


     do iq=1,nquad  ! loop over the quadrature points in each element
!!!!!!!!!!!!!!!!!!!!calculate shape function!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (nsd==2) then
	   if(nen_den.eq.3) then
		include 'sh2d3n.h'
	   elseif(nen_den.eq.4) then
		include 'sh2d4n.h'
	   end if
	elseif(nsd==3) then
	   if(nen_den.eq.4) then
		include 'sh3d4n.h'
	   elseif(nen_den.eq.8) then
		include 'sh3d8n.h'
	   end if
	end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     eft0 = abs(det)*wq(iq)   ! calculate weight at each quad point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     dr(:) = 0.0

     do inl=1,nen_den
	dr(1:nsd) = dr(1:nsd)+sh(1:nsd,inl)*d(inl)
     end do

     ph(0:nsd,1:nen_den) = sh(0:nsd,1:nen_den)*eft0

     do inl=1,nen_den
	node=ien_den(inl,ie)
	do isd=1,nsd
	   p_inter(node)=p_inter(node)+ph(isd,inl)*dr(isd)
	end do! residual
     end do

     end do ! end of quad loop

  end do ! end of ele loop
!write(*,*)'w_inter=',w_inter(:)
end subroutine blockgmres_Poisson










































