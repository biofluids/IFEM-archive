!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!residual for Laplace eqn!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  block_Laplace(xloc,I_fluid_center,p_inter,w_inter,ien_center)
  use global_constants
  use run_variables
  use fluid_variables,only:nsd,iquad
  use centermesh_variables
  implicit none

  integer ien_center(nen_center,ne_center)
  real(8) xloc(nsd,nn_center)
  real(8) I_fluid_center(nn_center)
  real(8) p_inter(nn_center)
  real(8) w_inter(nn_center)

  real(8) x(nsd,nen_center),d(nen_center)
  real(8) eft0,det,effd,effm,effc
  real(8) sh(0:nsd,nen_center),ph(0:nsd,nen_center)
  real(8) xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real(8) dr(nsd)
  real(8) hg
  
  integer inl,ie,isd,node,iq

  integer,parameter :: ndfpad=5,nsdpad=3,nenpad=8,nquadpad=8
  integer nquad
  real(8) sq(0:nsdpad,nenpad,nquadpad)
  real(8) xq(nsdpad,nquadpad),wq(nquadpad)

  if(nsd==3) then
     write(*,*)'no 3d in block laplace'
     stop
  end if

  if (nen_center==3) then
     call quad2d3n(iquad,nquad,xq,wq,nsdpad,nquadpad)
  end if
  do iq=1,nquad
     if(nen_center==3) then
        sq(0,1,iq) = xq(1,iq)
        sq(0,2,iq) = xq(2,iq)
        sq(0,3,iq) = 1 - xq(1,iq) - xq(2,iq)
     end if
  end do
!==============================================================
  do ie=1,ne_center     !loop over elements
     do inl=1,nen_center
	x(1:nsd,inl) = xloc(1:nsd,ien_center(inl,ie))
	d(inl) = I_fluid_center(ien_center(inl,ie))
     end do


     do iq=1,nquad  ! loop over the quadrature points in each element
!!!!!!!!!!!!!!!!!!!!calculate shape function!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (nsd==2) then
	   if(nen_center.eq.3) then
		include 'sh2d3n.h'
	   elseif(nen_center.eq.4) then
		include 'sh2d4n.h'
	   end if
	elseif(nsd==3) then
	   if(nen_center.eq.4) then
		include 'sh3d4n.h'
	   elseif(nen_center.eq.8) then
		include 'sh3d8n.h'
	   end if
	end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     eft0 = abs(det)*wq(iq)   ! calculate weight at each quad point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     dr(:) = 0.0

     do inl=1,nen_center
	dr(1:nsd) = dr(1:nsd)+sh(1:nsd,inl)*d(inl)
     end do

     ph(0:nsd,1:nen_center) = sh(0:nsd,1:nen_center)*eft0

     do inl=1,nen_center
	node=ien_center(inl,ie)
	do isd=1,nsd
	   p_inter(node)=p_inter(node)-ph(isd,inl)*dr(isd)
	end do! residual
	do isd=1,nsd
	   w_inter(node)=w_inter(node)+ph(isd,inl)*sh(isd,inl)
	end do! diag
     end do

     end do ! end of quad loop

  end do ! end of ele loop
!write(*,*)'w_inter=',w_inter(:)
end subroutine block_Laplace










































