!!!!!!!!get interface points normal using finite element!!!!!!!

subroutine get_inter_normal(xg,x_inter,ien,infdomain_inter,norm_fluid,norm_inter)

  use fluid_variables, only:nn,nsd,nen,ne
  use interface_variables

  real(8) xg(nsd,nn)
  real(8) x_inter(nsd,maxmatrix)
  integer ien(nen,ne)
  integer infdomain_inter(maxmatrix)
  real(8) norm_fluid(nsd,nn)
  real(8) norm_inter(nsd,maxmatrix)

  real(8) xg_ele(nsd,nen)
  integer node,inl,i,ie
  real(8) norm_ele(nsd,nen)
  real(8) sh(nen)
  real(8) xp(nsd)

  norm_inter(:,:) = 0.0
  do ie=1,ne
     do i=1,nn_inter
	if(infdomain_inter(i)==ie) then
	  do inl=1,nen
	     node=ien(inl,ie)
	     xg_ele(1:nsd,inl)=xg(1:nsd,node)
	     norm_ele(1:nsd,inl)=norm_fluid(1:nsd,node)
	  end do
	  xp(1:nsd)=x_inter(1:nsd,i)
	  call sh_exchange(xp,xg_ele,nsd,nen,sh)
	  do inl=1,nen
	     norm_inter(1:nsd,i)=norm_inter(1:nsd,i)+sh(inl)*norm_ele(1:nsd,inl)
	  end do
	end if
     end do
  end do

end subroutine get_inter_normal















