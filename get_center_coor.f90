!========================================
! get element center coordinates
!========================================

subroutine get_center_coor(xloc,x_center,ien,ien_center)

  use fluid_variables,only:nsd,nn,ne,nen
  use centermesh_variables, only:flag_center,nn_center,ne_center,nen_center

  real(8) xloc(nsd,nn)
  real(8) x_center(nsd,ne)
  integer ien(nen,ne)
  integer ien_center(nen_center,ne_center)

  integer i,j,ie,inl,isd,node
  real(8) temp(nsd,nen)
!  real(8) temp
  real(8) sh(nen)



!  if(flag_center==0) then
    x_center(:,:) = 0.0

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
!       do isd=1,nsd 
!	  temp = 0.0
!	  do inl=1,nen
!	     node=ien(inl,ie)
!	     temp=temp+xloc(isd,node)
!	  end do
!	  x_center(isd,ie)=temp/4.0
!       end do
    end do
   if(flag_center==0)then
    open(119,file='xyz_center.dat',status='unknown')
    do ie=1,ne
    write(119,999)x_center(1,ie),x_center(2,ie)
    end do
    close(119)
    999 format(f14.10,f14.10)
    write(*,*)'generate center coordiantes for matlab to generate mesh'
    stop
  elseif (flag_center==1) then
    open(118,file='ien_center.dat',status='old')
	do ie=1,ne_center
	   read(118,110)ien_center(1,ie),ien_center(2,ie),ien_center(3,ie)
	end do
    close(118)
110 format(I8,I8,I8)
  end if
!write(*,*)ien_center(:,:)
end subroutine get_center_coor
     
