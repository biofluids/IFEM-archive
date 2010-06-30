!=======================================================
!    residual for norm eqn
!=======================================================

subroutine block_norm(xloc,w,p_norm,ien,hk,I_fluid)

  use global_constants
  use run_variables
  use fluid_variables
  implicit none

  integer ien(nen,ne)
  real(8) xloc(nsd,nn)
  real(8) p_norm(nsd,nn)
  real(8) w(nsd,nn)
  real(8) I_fluid(nn)
  real(8) hk(ne)

  real(8) x(nsd,nen),d_I(nen)
  real(8) eft0,det,effd,effm,effc
  real(8) sh(0:nsd,nen),ph(0:nsd,nen)
  real(8) xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real(8) hg
  integer inl,ie,isd,node,iq
  real(8) dr_I(nsd)
  real(8) norm_I
!===========================================================
  do ie=1,ne     !loop over elements
     do inl=1,nen
        x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
        d_I(inl) = I_fluid(ien(inl,ie))
     end do

     hg = hk(ie)

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

	dr_I(:) = 0.0
	do inl=1,nen
	   dr_I(1:nsd) = dr_I(1:nsd)+sh(0,inl)*d_I(inl)
	end do
	ph(0:nsd,1:nen) = sh(0:nsd,1:nen)*eft0

	do inl=1,nen
	   node=ien(inl,ie)
	   do isd=1,nsd
		p_norm(isd,node)=p_norm(isd,node)+ph(isd,inl)*dr_I(isd)
		w(isd,node)=w(isd,node)+ph(0,inl)*sh(0,inl)
	   end do

	end do

     end do ! end of quad loop

  end do ! end of element loop
!write(*,*)'w_norm=',w(:,:)

end subroutine block_norm


























