!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!residual for Laplace eqn!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine block_curv(flag_domain,xloc,sur_fluid,curv_nn,p,w,ien,ne_local,ien_local,I_fluid)
  use global_constants
  use run_variables
  use fluid_variables
!  use centermesh_variables
!  use denmesh_variables, only:nn_den,ne_den,nen_den
  implicit none

  integer flag_domain(ne),ien(nen,ne)
  real(8) xloc(nsd,nn),sur_fluid(nsd,nn),curv_nn(nn)
  real(8) p(nn),w(nn),I_fluid(nn)


  real(8) x(nsd,nen),d(nen)
  real(8) eft0,det,effd,effm,effc
  real(8) sh(0:nsd,nen),ph(0:nsd,nen)
  real(8) xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real(8) dr(nsd)
  real(8) hg
  real(8) dI(nsd,nen),drI(nsd),II(nen)
  
  integer inl,ie,isd,node,iq,ine

  integer ne_local,ien_local(ne_local),ie_local


!==============================================================
!  do ie=1,ne_den     !loop over elements
!  do ine=1,ne_den_domain
  do ie_local=1,ne_local
     ie=ien_local(ie_local)
     if(flag_domain(ie)==0) goto 1111

     do inl=1,nen
	x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
	d(inl) = curv_nn(ien(inl,ie))
	dI(1:nsd,inl)=sur_fluid(:,ien(inl,ie))
	II(inl)=I_fluid(ien(inl,ie))
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
     drI(:)=0.0
     do inl=1,nen
	dr(1:nsd) = dr(1:nsd)+sh(1:nsd,inl)*d(inl)
     end do
     do inl=1,nen
	drI(1:nsd)=drI(1:nsd)+sh(1:nsd,inl)*II(inl)!*dI(1:nsd,inl)
     end do

     ph(0:nsd,1:nen) = sh(0:nsd,1:nen)*eft0

     do inl=1,nen
	node=ien(inl,ie)
	do isd=1,nsd
	   p(node)=p(node)-ph(0,inl)*dr(isd)*drI(isd)
	end do! residual

	do isd=1,nsd
	   w(node)=w(node)+ph(0,inl)*sh(isd,inl)*drI(isd)
	end do! diag
     end do

     end do ! end of quad loop

1111 continue
  end do ! end of ele loop
end subroutine block_curv










































