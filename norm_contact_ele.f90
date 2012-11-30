

subroutine norm_contact_ele(nn_con_ele,con_ele,norm_con_ele,ien,norm_node)

  use fluid_variables,only:nsd,nn,ne,nen
  use mpi_variables

  integer nn_con_ele,con_ele(nn_con_ele)
  real(8) norm_con_ele(nsd,nn_con_ele),norm_node(nsd,nn)
  integer ien(nen,ne)

  integer i,j,icount,jcount,isd,inl,node,ie
  real(8) temp(nsd),err_p
  real(8) modmod
  err_p=1.0e-4


  do i=1,nn_con_ele
     ie=con_ele(i)
     icount=0
     temp(:)=0.0
     do inl=1,nen
	node=ien(inl,ie)
if(nsd==2) then
	if( (abs(norm_node(1,node)).gt.err_p) .or. &
	    (abs(norm_node(2,node)).gt.err_p)) then
	  icount=icount+1
	  temp(1:nsd)=temp(1:nsd)+norm_node(1:nsd,node)
	end if

else

        if( (abs(norm_node(1,node)).gt.err_p) .or. &
            (abs(norm_node(2,node)).gt.err_p) .or. &
            (abs(norm_node(3,node)).gt.err_p)) then
          icount=icount+1
          temp(1:nsd)=temp(1:nsd)+norm_node(1:nsd,node)
        end if



end if

     end do
     temp(1:nsd)=temp(1:nsd)/real(icount)
     modmod=0
     do j=1,nsd
        modmod=modmod+temp(j)**2
     end do
     modmod=sqrt(modmod)

     temp(1:nsd)=temp(1:nsd)/modmod!/(sqrt(temp(1)**2+temp(2)**2))
     norm_con_ele(1:nsd,i)=-temp(1:nsd)
!    norm_con_ele(1,i)=1.0
!    norm_con_ele(2:nsd,i)=1.0
  end do

end subroutine norm_contact_ele
  
