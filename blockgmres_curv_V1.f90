!=======================================================
!    residual for norm eqn
!=======================================================

subroutine blockgmres_curv(xloc,norm_fluid_sub,p_norm,ien,hk)

  use global_constants
  use run_variables
  use fluid_variables
  implicit none

  integer ien(nen,ne)
  real(8) xloc(nsd,nn)
  real(8) norm_fluid_sub(nn)
  real(8) p_norm(nn)
  real(8) hk(ne)

  real(8) x(nsd,nen),d_I(nen)
  real(8) eft0,det,effd,effm,effc
  real(8) sh(0:nsd,nen),ph(0:nsd,nen)
  real(8) xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real(8) hg
  integer inl,ie,isd,node,iq
  real(8) ds_I
!===========================================================
  do ie=1,ne     !loop over elements
     do inl=1,nen
        x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
!        d_I(inl) = I_fluid(ien(inl,ie))
	d_I(inl) = norm_fluid_sub(ien(inl,ie))
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

	ds_I = 0.0
	do inl = 1,nen
	   ds_I = ds_I+sh(0,inl)*d_I(inl)
	end do

	ph(0:nsd,1:nen) = sh(0:nsd,1:nen)*eft0

	do inl=1,nen
	   node=ien(inl,ie)
		p_norm(node)=p_norm(node)+ph(0,inl)*ds_I

	end do

     end do ! end of quad loop

  end do ! end of element loop

end subroutine blockgmres_curv


























