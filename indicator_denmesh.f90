!===============================================================!
! get the indicator for denmesh                 !
!===============================================================!

subroutine indicator_denmesh(I_fluid_den,x_den,ien_den,nn_den,ne_den,nen_den,I_fluid_center)

  use mpi_variables
!  use denmesh_variables  !nn_den,ne_den,nen_den,nbc_den
  use fluid_variables, only:nsd
  use allocate_variables, only:ne_inter,inter_ele, &
				den_domain,ne_den_domain
!inter_ele_den(1:ne_inter_den),den_domain(1:ne_den_domain),inter_ele_den(1:ne_inter_den)

  integer nn_den,ne_den,nen_den,nbc_den
  real(8) I_fluid_den(nn_den),I_fluid_center(ne_den),I_fluid_temp(nn_den)
  real(8) x_den(nsd,nn_den)
  integer ien_den(nen_den,ne_den)
!  integer bcnode_den(nbc_den,2)
!  integer its
  integer  set_flag

  real(8) p(nn_den),w(nn_den) !residual and preconditioner
  real(8) dg(nn_den)          !solution
  integer offdomain(ne_den),ne_offdomain  !elements that are not in the calculation domain
!  integer flag_offdomain(ne_den)

  integer fake_domain(ne_den),fake_ne
  integer i,j,icount,jcount,ie,node,inl
  I_fluid_temp=I_fluid_den
  p(:)=0.0
  w(:)=0.0
    ne_offdomain=0
!    flag_offdomain(:)=0
    fake_ne=0
    do i=1,ne_den
       set_flag=0
       do j=1,nen_den
          if( (I_fluid_den(ien_den(j,i)).gt.1.0e-6) .and. (I_fluid_den(ien_den(j,i)).lt.1.0-1.0e-6)) then
             set_flag=1
          end if
        end do

!       if((I_fluid_center(i).gt.1.0e-6).and.(I_fluid_center(i).lt.(1.0-1.0e-6))) then
        if(set_flag==1) then
!         flag_offdomain(i)=1
	 fake_ne=fake_ne+1
         fake_domain(fake_ne)=i
	else
	 ne_offdomain=ne_offdomain+1
	 offdomain(ne_offdomain)=i
	end if
    end do



!    do i=1,ne_den_domain
!	flag_offdomain(den_domain(i))=1
!    end do
!    do i=1,ne_den
!       if(flag_offdomain(i)==0) then
!	  ne_offdomain=ne_offdomain+1
!	  offdomain(ne_offdomain)=i
!       end if
!    end do  !find out offdomian
!    fake_ne=ne_den_domain
!    fake_domain(1:fake_ne)=den_domain(1:ne_den_domain)

	   
	
!==============================================================
!  set boundary condition
!    do i=1,nbc_den
!	I_fluid_den(bcnode_den(i,1))=real(bcnode_den(i,2))
!    end do
    do i=1,ne_inter
	ie=inter_ele(i)
	do inl=1,nen_den
	   node=ien_den(inl,ie)
	   I_fluid_den(node)=0.5
	end do
    end do           
!==============================================================
    call block_Laplace(fake_ne,fake_domain,x_den,I_fluid_den,p,w,ien_den,nn_den,ne_den,nen_den)

!==============================================================
!  set residual at bc and offdomain to be 0
    do i=1,ne_offdomain
       ie=offdomain(i)
       do inl=1,nen_den
	  node=ien_den(inl,ie)
	  p(node)=0.0
       end do
    end do
!    do i=1,nbc_den
!        p(bcnode_den(i,1))=0.0
!    end do
    do i=1,ne_inter
        ie=inter_ele(i)
        do inl=1,nen_den
           node=ien_den(inl,ie)
           p(node)=0.0
        end do
    end do
!==============================================================
    dg(:)=0.0
    call gmres_Laplace(fake_ne,fake_domain,x_den,I_fluid_den,w,p,dg,ien_den,&
			offdomain,ne_offdomain,nn_den,ne_den,nen_den)


   I_fluid_den(:)=I_fluid_den(:)+dg(:)
   do i=1,nn_den
      if(I_fluid_den(i).gt.0.5001)then
        I_fluid_den(i)=1.0
      else  if(I_fluid_den(i).lt.0.4999) then
        I_fluid_den(i)=0.0
      end if
   end do
   do i=1,ne_den
      I_fluid_center(i)=0.0
      do inl=1,nen_den
         I_fluid_center(i)=I_fluid_center(i)+1.0/nen_den*I_fluid_den(ien_den(inl,i))
      end do
      if(I_fluid_center(i).ge.0.5) I_fluid_center(i)=1.0
      if(I_fluid_center(i).lt.0.5) I_fluid_center(i)=0.0
    end do
   
    do icount=1,ne_den_domain
        ie=den_domain(icount)
	I_fluid_center(ie)=0.0
	do inl=1,nen_den
	   I_fluid_center(ie)=I_fluid_center(ie)+1.0/real(nen_den)*I_fluid_temp(ien_den(inl,ie))
	end do
    end do
end subroutine indicator_denmesh






