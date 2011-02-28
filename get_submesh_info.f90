!========================================
! get element center coordinates
!========================================

subroutine get_submesh_info(xloc,x_center,ien,ien_den,x_den,bcnode_den)

  use fluid_variables,only:nsd,nn,ne,nen
  use denmesh_variables, only:flag_den,nn_den,ne_den,nen_den,nbc_den

  real(8) xloc(nsd,nn)
  real(8) x_center(nsd,ne)
  integer ien(nen,ne)
  integer ien_den(nen_den,ne_den)
  real(8) x_den(nsd,nn_den)
  integer bcnode_den(nbc_den)

  integer i,j,ie,inl,isd,node
  real(8) temp(nsd,nen)
!  real(8) temp
  real(8) sh(nen)



!  if(flag_center==0) then
    x_center(:,:) = 0.0
!get center coordinates
    do ie = 1,ne
	if((nsd==2).and.(nen==4))then
	  sh(:)=0.25
	  do inl=1,nen
	     node=ien(inl,ie)
	     temp(1:nsd,inl)=xloc(1:nsd,node)
	  end do
	  do inl=1,nen
	     x_center(1:nsd,ie)=x_center(1:nsd,ie)+sh(inl)*temp(1:nsd,inl)
	  end do
	end if
    end do
if(nen_den.ne.4) then
  write(*,*)'dense_mesh should be rectangular right now'
  stop
end if
    open(118,file='ien_den.dat',status='old')
	do ie=1,ne_den
!	   read(118,110)ien_den(1,ie),ien_den(2,ie),ien_den(3,ie)
!	  read(118,110)ien_den(1,ie),ien_den(2,ie),ien_den(3,ie),ien_den(4,ie)
        read(118,110)ien_den(1:nen_den,ie)
	end do
    close(118)
    open(119,file='xyz_den.dat',status='old')
	do i=1,nn_den
!	   read(119,111)x_center_den(1,i), x_center_den(2,i)
	  read(119,111)x_den(1:nsd,i)
	end do
    close(119)
    open(120,file='bcnode_den.dat',status='old')
	do i=1,nbc_den
	   read(120,112)bcnode_den(i)
	end do
    close(12)
110 format(I8,I8,I8,I8)
111 format(f14.10,f14.10)
112 format(I8)
!  end if
!write(*,*)ien_center(:,:)
    if(nsd==3) then
	write(*,*)'no 3d for centermesh'
	stop
     end if
    
end subroutine get_submesh_info
     
