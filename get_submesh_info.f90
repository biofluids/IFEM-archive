!========================================
! get element center coordinates
!========================================

subroutine get_submesh_info(xloc,x_center,ien,bcnode)

  use fluid_variables,only:nsd,nn,ne,nen
!  use denmesh_variables, only:nn_den,ne_den,nen_den,nbc_den
  use interface_variables, only:nbc
  real(8) xloc(nsd,nn)
  real(8) x_center(nsd,ne)
  integer ien(nen,ne)
!  integer ien_den(nen_den,ne_den)
!  real(8) x_den(nsd,nn_den)
  integer bcnode(nbc,2)

  integer i,j,ie,inl,isd,node
  real(8) temp(nsd,nen)
  real(8) sh(nen)



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


!    open(120,file='bcnode.dat',status='old')
!	do i=1,nbc
!	   read(120,112)bcnode(i,1),bcnode(i,2)
!	end do
!    close(120)
112 format(2I8)
    
end subroutine get_submesh_info
     
