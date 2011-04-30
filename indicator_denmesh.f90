!===============================================================!
! get the indicator for denmesh                 !
!===============================================================!

subroutine indicator_denmesh(I_fluid_den,x_den,ien_den,bcnode_den,its)

  use mpi_variables
  use denmesh_variables  !nn_den,ne_den,nen_den,nbc_den
  use fluid_variables, only:nsd
  use allocate_variables, only:ne_inter_den,inter_ele_den, &
				den_domain,ne_den_domain
!inter_ele_den(1:ne_inter_den),den_domain(1:ne_den_domain),inter_ele_den(1:ne_inter_den)

  real(8) I_fluid_den(nn_den)
  real(8) x_den(nsd,nn_den)
  integer ien_den(nen_den,ne_den)
  integer bcnode_den(nbc_den,2)
  integer its

  real(8) p(nn_den),w(nn_den) !residual and preconditioner
  real(8) dg(nn_den)          !solution
  integer offdomain(ne_den),ne_offdomain  !elements that are not in the calculation domain
  integer flag_offdomain(ne_den)

  integer fake_domain(ne_den),fake_ne
  integer i,j,icount,jcount,ie,node,inl

  p(:)=0.0
  w(:)=0.0
  if(its==1) then
    I_fluid_den(:)=0.0
    ne_offdomain=0
    offdomain(:)=0
    fake_ne=ne_den
    do i=1,ne_den
       fake_domain(i)=i
    end do            !it's actually the entire domain
  else
    ne_offdomain=0
    flag_offdomain(:)=0
    do i=1,ne_den_domain
	flag_offdomain(den_domain(i))=1
    end do
    do i=1,ne_den
       if(flag_offdomain(i)==0) then
	  ne_offdomain=ne_offdomain+1
	  offdomain(ne_offdomain)=i
       end if
    end do  !find out offdomian
    fake_ne=ne_den_domain
    fake_domain(1:fake_ne)=den_domain(1:ne_den_domain)

  end if
	   
	
!==============================================================
!  set boundary condition
    do i=1,nbc_den
	I_fluid_den(bcnode_den(i,1))=real(bcnode_den(i,2))
    end do
    do i=1,ne_inter_den
	ie=inter_ele_den(i)
	do inl=1,nen_den
	   node=ien_den(inl,ie)
	   I_fluid_den(node)=1.0
	end do
    end do           
!==============================================================
    call block_Laplace(fake_ne,fake_domain,x_den,I_fluid_den,p,w,ien_den)

!==============================================================
!  set residual at bc and offdomain to be 0
    do i=1,ne_offdomain
       ie=offdomain(i)
       do inl=1,nen_den
	  node=ien_den(inl,ie)
	  p(node)=0.0
       end do
    end do
    do i=1,nbc_den
        p(bcnode_den(i,1))=0.0
    end do
    do i=1,ne_inter_den
        ie=inter_ele_den(i)
        do inl=1,nen_den
           node=ien_den(inl,ie)
           p(node)=0.0
        end do
    end do
!==============================================================
    dg(:)=0.0
    call gmres_Laplace(fake_ne,fake_domain,x_den,I_fluid_den,w,p,dg,ien_den,bcnode_den,&
			offdomain,ne_offdomain)


   I_fluid_den(:)=I_fluid_den(:)+dg(:)
if(myid==0) then
!write(*,*)'I_fluid_den=',I_fluid_den(:)
end if



end subroutine indicator_denmesh






