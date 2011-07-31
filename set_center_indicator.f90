!===============================================
!setcenter indicator
!==================================================

subroutine set_center_indicator(I_fluid_den,I_fluid_center,ien_den,its,infdomain_den,flag_den,&
				ien,I_fluid)
  use mpi_variables
  use interface_variables
  use denmesh_variables
  use fluid_variables, only:ne,nen,nn
  use allocate_variables, only:den_domain,center_domain, &
				ne_den_domain,nn_center_domain, &
				ne_inter_den,inter_ele_den,&
				ne_inter,inter_ele, &
				ne_regen_ele,regen_ele

  real(8) I_fluid_den(nn_den)
  real(8) I_fluid_center(ne)
  integer ien_den(nen_den,ne_den)
  integer its
  integer infdomain_den(ne)
  real(8) I_fluid(nn)
  integer ien(nen,ne)
  

  real(8) eps
  integer ie,je,icount,inl,node,flag
  integer flag_den(ne_den)

  eps=0.8e-4
  if (its==1) then
    flag_den(:)=1 !set initial flag of den element to be 1 inner element
    do ie=1,ne_den
	flag=0
	do inl=1,nen_den
	   node=ien_den(inl,ie)
	   if(abs(I_fluid_den(node)-1.0).gt.eps) then
	     flag=1
	     goto 101
	   end if
	end do
101 continue
	if(flag==1) then
	   flag_den(ie)=-1
	end if   ! find out outer element and set flag to be -1
    end do

    do ie=1,ne_inter_den
	flag_den(inter_ele_den(ie))=0
    end do      ! set inter element flag to be 0
  else
!    flag_den(:)=1
    do ie=1,ne_den_domain
	flag_den(den_domain(ie))=1
	flag=0
	do inl=1,nen_den
	   node=ien_den(inl,den_domain(ie))
	   if(abs(I_fluid_den(node)-1.0).gt.eps) then
	     flag=1
	     goto 102
	   end if
	end do
102 continue
	if(flag==1) then
	  flag_den(den_domain(ie))=-1
	end if
    end do
    do ie=1,ne_inter_den
	flag_den(inter_ele_den(ie))=1
    end do
  end if
!==================================================
  if (its==1) then
    do ie=1,ne
	if(flag_den(infdomain_den(ie))==1) then
	  I_fluid_center(ie)=1.0
	else if(flag_den(infdomain_den(ie))==0) then
	  I_fluid_center(ie)=0.5
	else
	  I_fluid_center(ie)=0.0
	end if
    end do
  else
       
!  do icount=1,nn_center_domain
!     ie=center_domain(icount)
!     if(flag_den(infdomain_den(ie))==1) then
!	I_fluid_center(ie)=1.0
!     else if(flag_den(infdomain_den(ie))==0) then
!	I_fluid_center(ie)=0.5
!     else
!	I_fluid_center(ie)=0.0
!     end if
!  end do   !set up initial guess for center point
!
!  do icount=1,ne_den_domain
!     ie=den_domain(icount)
!     I_fluid_center(ie)=0.0
!     do inl=1,nen
!	I_fluid_center(ie)=I_fluid_center(ie)+0.25*I_fluid(ien(inl,ie))
!     end do
!  end do   !set up inter center indicator using last time step I_fluid

  end if !======================================

  if (its==1) then
    do ie=1,ne_den
	if(flag_den(ie)==1) then
	   I_fluid_den(ien_den(1:nen_den,ie))=1.0
	else! if(flag_den(ie)==-1) then
	   I_fluid_den(ien_den(1:nen_den,ie))=0.0
	end if
    end do
  else
    do icount=1,ne_den_domain
	ie=den_domain(icount)
	if(flag_den(ie)==1) then
	   I_fluid_den(ien_den(1:nen_den,ie))=1.0
	else! if(flag_den(ie)==-1) then
	   I_fluid_den(ien_den(1:nen_den,ie))=0.0
	end if
    end do
  end if     ! update I_fluid_den

end subroutine set_center_indicator













