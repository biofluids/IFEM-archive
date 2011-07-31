

subroutine indicator_denmesh_after(x_inter,I_fluid_den,x_den,ien_den,bcnode_den,norm_inter,arc_inter,vol_nn)

  use mpi_variables
  use interface_variables
  use denmesh_variables !nn_den,ne_den,nen_den,nbc_den
  use fluid_variables, only:nsd
  use allocate_variables, only:den_domain,ne_den_domain

  real(8) I_fluid_den(nn_den)
  real(8) x_den(nsd,nn_den)
  real(8) x_inter(nsd,maxmatrix)
  integer ien_den(nen_den,ne_den)
  integer bcnode_den(nbc_den,2)
  real(8) norm_inter(nsd,maxmatrix),arc_inter(maxmatrix)
  real(8) vol_nn(nn_den)

  real(8) p(nn_den),w(nn_den)
  real(8) dg(nn_den)
  integer offdomain(ne_den),ne_offdomain
  integer flag_offdomain(ne_den)

  integer fake_domain(ne_den),fake_ne
  integer i,j,icount,jcount,ie,node,inl

  real(8) norm_fluid(nsd,nn_den)

  call norm_distribution(x_den,x_inter,vol_nn,arc_inter,norm_inter,norm_fluid)
if(myid==0) write(*,*)'finish norm distribution'
  p(:)=0.0
  w(:)=0.0
  ne_offdomain=0
  flag_offdomain(:)=0
  do i=1,ne_den_domain
      flag_offdomain(den_domain(i))=1
  end do
  do i=1,ne_den
     if(flag_offdomain(i)==0) then
	do inl=1,nen_den
	   if(I_fluid_den(ien_den(inl,i)).ge.0.5) then
	      I_fluid_den(ien_den(inl,i))=1.0
	   else
	      I_fluid_den(ien_den(inl,i))=0.0
	   end if
	end do
        ne_offdomain=ne_offdomain+1
        offdomain(ne_offdomain)=i
     end if
  end do  !find out offdomian
  fake_ne=ne_den_domain
  fake_domain(1:fake_ne)=den_domain(1:ne_den_domain)

  call block_Poisson(fake_ne,fake_domain,x_den,I_fluid_den,p,w,ien_den,norm_fluid)
if(myid==0) write(*,*)'finish block poisson'
  do i=1,ne_offdomain
     ie=offdomain(i)
     do inl=1,nen_den
        node=ien_den(inl,ie)
        p(node)=0.0
     end do
  end do

  dg(:)=0.0
  call gmres_Poisson(fake_ne,fake_domain,x_den,I_fluid_den,w,p,dg,ien_den,bcnode_den,&
                      offdomain,ne_offdomain)

  I_fluid_den(:)=I_fluid_den(:)+dg(:)

end subroutine indicator_denmesh_after



