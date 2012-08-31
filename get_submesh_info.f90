!========================================
! get element center coordinates
!========================================

subroutine get_submesh_info(xloc,x_center,ien,vol_nn,hg)

  use fluid_variables,only:nsd,nn,ne,nen
  real(8) xloc(nsd,nn)
  real(8) x_center(nsd,ne)
  integer ien(nen,ne)

  integer i,j,ie,inl,isd,node
  real(8) temp(nsd,nen)
  real(8) sh(nen)

  real(8) vol_nn(nn),hg(ne)


    x_center(:,:) = 0.0
!get center coordinates
    do ie = 1,ne
	if((nsd==2).and.(nen==4))then
	  sh(:)=0.25        !2d4n
	else if(nsd==3.and.nen==8) then
	  sh(:)=0.125       !3d8n
	end if
	  do inl=1,nen
	     node=ien(inl,ie)
	     temp(1:nsd,inl)=xloc(1:nsd,node)
	  end do
	  do inl=1,nen
	     x_center(1:nsd,ie)=x_center(1:nsd,ie)+sh(inl)*temp(1:nsd,inl)
	  end do
    end do

   vol_nn(:)=0.0
   do ie=1,ne
      vol_nn(ien(1:nen,ie))=vol_nn(ien(1:nen,ie))+hg(ie)**nsd/real(nen)
   end do
end subroutine get_submesh_info
     
