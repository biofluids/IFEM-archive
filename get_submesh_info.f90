!========================================
! get element center coordinates
!========================================

subroutine get_submesh_info(xloc,x_center,ien,ien_den,x_den,bcnode_den)

  use fluid_variables,only:nsd,nn,ne,nen
  use denmesh_variables, only:nn_den,ne_den,nen_den,nbc_den

  real(8) xloc(nsd,nn)
  real(8) x_center(nsd,ne)
  integer ien(nen,ne)
  integer ien_den(nen_den,ne_den)
  real(8) x_den(nsd,nn_den)
  integer bcnode_den(nbc_den,2)

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




    open(118,file='ien_den.dat',status='old')
      if(ne_den==3) then
	do ie=1,ne_den
        read(118,110)ien_den(1:nen_den,ie)
	end do
      else if(nen_den==4) then
	do ie=1,ne_den
	read(118,100)ien_den(1:nen_den,ie)
	end do
      end if
    close(118)
110 format(3I8)
100 format(4I8)


    open(119,file='xyz_den.dat',status='old')
      if(nsd==2) then
	do i=1,nn_den
	  read(119,111)x_den(1:nsd,i)
	end do
      else if(nsd==3) then
	do i=1,nn_den
	  read(119,101)x_den(1:nsd,i)
	end do
      end if
    close(119)
111 format(2f14.10)
101 format(3f14.10)


    open(120,file='bcnode_den.dat',status='old')
	do i=1,nbc_den
	   read(120,112)bcnode_den(i,1),bcnode_den(i,2)
	end do
    close(120)
112 format(2I8)
    
end subroutine get_submesh_info
     
