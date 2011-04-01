!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!residual for Laplace eqn!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  block_Laplace(xloc,I_fluid,p_inter,w_inter,ien,lp_el,count_el,source)
  use global_constants
  use run_variables
  use fluid_variables
  implicit none

  integer ien(nen,ne)
  real(8) xloc(nsd,nn)
  real(8) I_fluid(nn)
  real(8) p_inter(nn)
  real(8) w_inter(nn)
!----------------------------
  integer count_el
  integer lp_el(count_el)
  real(8) source(nsd,nn)
!----------------------------
  real(8) x(nsd,nen),d(nen)
  real(8) eft0,det,effd,effm,effc
  real(8) sh(0:nsd,nen),ph(0:nsd,nen)
  real(8) xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)
  
  real(8) dr(nsd)
  
  integer inl,ie,isd,node,iq
!---------------------------------
  integer je
  real(8) local_source(nsd,nen)
  real(8) dg_l(nsd)
  real(8) divg
!==============================================================
  do je=1,count_el     !loop over elements
	ie=lp_el(je)
!   do ie=1,ne
     do inl=1,nen
	x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
	d(inl) = I_fluid(ien(inl,ie))
	local_source(1:nsd,inl)=source(1:nsd,ien(inl,ie))
     end do


     do iq=1,nquad  ! loop over the quadrature points in each element
!!!!!!!!!!!!!!!!!!!!calculate shape function!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (nsd==2) then
	   if(nen.eq.3) then
		include 'sh2d3n.h'
	   elseif(nen.eq.4) then
		include 'sh2d4n.h'
	   end if
	elseif(nsd==3) then
	   if(nen.eq.4) then
		include 'sh3d4n.h'
	   elseif(nen.eq.8) then
		include 'sh3d8n.h'
	   end if
	end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     eft0 = abs(det)*wq(iq)   ! calculate weight at each quad point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     dr(:) = 0.0
     dg_l(:)=0.0

     do inl=1,nen
	dr(1:nsd) = dr(1:nsd)+sh(1:nsd,inl)*d(inl)
	dg_l(1:nsd) = dg_l(1:nsd)+sh(1:nsd,inl)*local_source(1:nsd,inl) ! dG/dx
     end do
	divg=sum(dg_l) ! div(G)

     ph(0:nsd,1:nen) = sh(0:nsd,1:nen)*eft0

     do inl=1,nen
	node=ien(inl,ie)
	do isd=1,nsd
	   p_inter(node)=p_inter(node)-ph(isd,inl)*dr(isd)+ph(0,inl)*divg
	end do! residual
	do isd=1,nsd
	   w_inter(node)=w_inter(node)+ph(isd,inl)*sh(isd,inl)
	end do! diag
     end do

     end do ! end of quad loop

  end do ! end of ele loop
!write(*,*)'w_inter=',w_inter(:)
end subroutine block_Laplace










































