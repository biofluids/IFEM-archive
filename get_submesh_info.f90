!========================================
! get element center coordinates
!========================================

subroutine get_submesh_info(xloc,x_center,ien,ne_spbc,spbcele,hg)

  use fluid_variables,only:nsd,nn,ne,nen
  use interface_variables, only:nn_center,c_w ! weight for the center points
  use mpi_variables
  real(8) xloc(nsd,nn)
  real(8) x_center(nsd,nn_center)
  integer ien(nen,ne)
  integer ne_spbc,spbcele(ne_spbc)
  real(8) hg(ne)

  integer i,j,k,ie,inl,isd,node
  real(8) temp(nsd,nen)
  real(8) sh(nen)
  integer flag,ccount

  integer, parameter :: nn_sub=2
  real(8) shap(nen,nn_sub**3),x_loc_can(nsd,nn_sub**3)

  if(allocated(c_w)) deallocate(c_w)
  allocate(c_w(nn_center))

    x_center(:,:) = 0.0
    ccount=0
!get center coordinates
    do ie = 1,ne
	flag=1
	do i=1,ne_spbc
	   if(spbcele(i)==ie) flag=0
	end do
	if(flag==1) then
	  if((nsd==2).and.(nen==4))then
	    sh(:)=0.25        !2d4n
	  else if(nsd==3.and.nen==8) then
	    sh(:)=0.125       !3d8n
	  end if
	  ccount=ccount+1
	  do inl=1,nen
	     node=ien(inl,ie)
	     temp(1:nsd,inl)=xloc(1:nsd,node)
	  end do
	  do inl=1,nen
	     x_center(1:nsd,ccount)=x_center(1:nsd,ccount)+sh(inl)*temp(1:nsd,inl)
	  end do
	  c_w(ccount)=hg(ie)
	end if
    end do

    do i=1,nn_sub
       do j=1,nn_sub
          do k=1,nn_sub
             node=nn_sub*(nn_sub*(i-1)+j-1)+k
             x_loc_can(1,nn_sub*(nn_sub*(i-1)+j-1)+k)=2.0/real(nn_sub)*i-1.0-1.0/real(nn_sub)
             x_loc_can(2,nn_sub*(nn_sub*(i-1)+j-1)+k)=2.0/real(nn_sub)*j-1.0-1.0/real(nn_sub)
             x_loc_can(3,nn_sub*(nn_sub*(i-1)+j-1)+k)=2.0/real(nn_sub)*k-1.0-1.0/real(nn_sub)
             shap(1,node)=0.125*(1-x_loc_can(1,node))*(1-x_loc_can(2,node))*(1-x_loc_can(3,node))
             shap(2,node)=0.125*(1+x_loc_can(1,node))*(1-x_loc_can(2,node))*(1-x_loc_can(3,node))
             shap(3,node)=0.125*(1+x_loc_can(1,node))*(1+x_loc_can(2,node))*(1-x_loc_can(3,node))
             shap(4,node)=0.125*(1-x_loc_can(1,node))*(1+x_loc_can(2,node))*(1-x_loc_can(3,node))
             shap(5,node)=0.125*(1-x_loc_can(1,node))*(1-x_loc_can(2,node))*(1+x_loc_can(3,node))
             shap(6,node)=0.125*(1+x_loc_can(1,node))*(1-x_loc_can(2,node))*(1+x_loc_can(3,node))
             shap(7,node)=0.125*(1+x_loc_can(1,node))*(1+x_loc_can(2,node))*(1+x_loc_can(3,node))
             shap(8,node)=0.125*(1-x_loc_can(1,node))*(1+x_loc_can(2,node))*(1+x_loc_can(3,node))

          end do
       end do
    end do


    do ie=1,ne_spbc
       do inl=1,nen
          node=ien(inl,spbcele(ie))
	  temp(1:nsd,inl)=xloc(1:nsd,node)
       end do
       do icount=1,nn_sub**nsd
	  ccount=ccount+1
          do inl=1,nen
	     x_center(1:nsd,ccount)=x_center(1:nsd,ccount)+shap(inl,icount)*temp(1:nsd,inl)
	  end do
	  c_w(ccount)=hg(spbcele(ie))/real(nn_sub**nsd)
	end do
    end do

    
end subroutine get_submesh_info
     
