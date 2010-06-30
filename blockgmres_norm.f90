!=======================================================
!    residual for norm eqn
!=======================================================

subroutine blockgmres_norm(xloc,norm_fluid_sub,p_norm,ien,hk)

  use global_constants
  use run_variables
  use fluid_variables
  implicit none

  integer ien(nen,ne)
  real(8) xloc(nsd,nn)
  real(8) norm_fluid_sub(nsd,nn)
  real(8) p_norm(nsd,nn)
  real(8) hk(ne)

  real(8) x(nsd,nen),d_I(nen),d_norm(nsd,nen)
  real(8) eft0,det,effd,effm,effc
  real(8) sh(0:nsd,nen),ph(0:nsd,nen)
  real(8) xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real(8) hg
  integer inl,ie,isd,node,iq
  real(8) dr_I(nsd)
  real(8) ds_norm(nsd)
  real(8) norm_I
!===========================================================
  do ie=1,ne     !loop over elements
     do inl=1,nen
        x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
!        d_I(inl) = I_fluid(ien(inl,ie))
	d_norm(1:nsd,inl) = norm_fluid_sub(1:nsd,ien(inl,ie))
     end do

     hg = hk(ie)
!write(*,*)'d_norm=',d_norm(:,:)
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

	eft0 = abs(det)*wq(iq)

!	dr_I(:) = 0.0
!	do inl=1,nen
!	   dr_I(1:nsd) = dr_I(1:nsd)+sh(1:nsd,inl)*d_I(inl)
!	end do

!	norm_I = 0.0
!	do isd = 1,nsd
!	   norm_I = norm_I+dr_I(isd)**2
!	end do
!	norm_I = sqrt(norm_I)
	ds_norm(:) = 0.0
	do inl = 1,nen
	   ds_norm(1:nsd) = ds_norm(1:nsd)+sh(0,inl)*d_norm(1:nsd,inl)
	end do

	ph(0:nsd,1:nen) = sh(0:nsd,1:nen)*eft0

	do inl=1,nen
	   node=ien(inl,ie)
	   do isd=1,nsd
		p_norm(isd,node)=p_norm(isd,node)+ph(0,inl)*ds_norm(isd)
	   end do

	end do

     end do ! end of quad loop

  end do ! end of element loop

end subroutine blockgmres_norm


























