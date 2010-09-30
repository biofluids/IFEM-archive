!=================================
!  residual for curv equation
!=================================

subroutine block_curv(xloc,w_curv,p_curv,ien,hk,I_fluid)

  use global_constants
  use run_variables
  use fluid_variables
  implicit none

  integer ien(nen,ne)
  real(8) xloc(nsd,nn)
  real(8) p_curv(nn)
  real(8) w_curv(nn)
  real(8) I_fluid(nn)
  real(8) hk(ne)

  real(8) x(nsd,nen),d_I(nen)
  real(8) eft0,det,effd,effm,effc
  real(8) sh(0:nsd,nen),ph(0:nsd,nen)
  real(8) xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real(8) hg
  integer inl,ie,isd,node,iq
  real(8) dr_I(nsd)
  real(8) temp
  real(8) eps

!=========================================
  eps=0.0001
  do ie=1,ne    ! loop over elements
     do inl=1,nen
	x(1:nsd,inl)=xloc(1:nsd,ien(inl,ie))
	d_I(inl)=I_fluid(ien(inl,ie))
     end do

     hg=hk(ie)

     do iq=1,nquad  ! loop over the quadrature points
!===============calculate shape function==================
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	dr_I(:)=0.0

	do inl=1,nen
	   dr_I(1:nsd)=dr_I(1:nsd)+sh(1:nsd,inl)*d_I(inl)
	end do

	temp=0.0
	do isd=1,nsd
	   temp=temp+dr_I(isd)**2
	end do
	temp=sqrt(temp)
	ph(0:nsd,1:nen)=sh(0:nsd,1:nen)*eft0

	if (temp.gt.1.0e-6) then

	do inl=1,nen
	   node=ien(inl,ie)
	   do isd=1,nsd
		p_curv(node)=p_curv(node)+ph(isd,inl)*dr_I(isd)/temp
	   end do
	end do

	end if

	do inl=1,nen
	   node=ien(inl,ie)
	   w_curv(node)=w_curv(node)+ph(0,inl)*sh(0,inl)
	   do isd=1,nsd
	      w_curv(node)=w_curv(node)+eps*ph(isd,inl)*sh(isd,inl)
	   end do
	end do

   end do ! end of quad loop

  end do ! end of ele loop
!write(*,*)'p_curv=',p_curv(:)
!write(*,*)'w_curv=',w_curv(:)

end subroutine block_curv











































