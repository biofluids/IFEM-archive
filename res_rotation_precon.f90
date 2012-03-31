

subroutine res_rotation_precon(p,norm_node,nn_spbc,spbcnode)

  use fluid_variables, only:ndf,nn,nsd,ne
  use mpi_variables

  real(8) p(ndf,nn)
  real(8) norm_node(nsd,nn)
  integer nn_spbc,spbcnode(nn_spbc)

  real(8) norm(nsd),tang(nsd)
  real(8) temp(nsd)
  integer i,j,node

  do i=1,nn_spbc
     node=spbcnode(i)
     norm(1:nsd)=norm_node(1:nsd,node)
     tang(1)=-norm(2)
     tang(2)=norm(1)
     temp(1)=norm(1)**2*p(1,node)+norm(2)**2*p(2,node)
     temp(2)=tang(1)**2*p(1,node)+tang(2)**2*p(2,node)
     p(1,node)=temp(1)
     p(2,node)=temp(2)
  end do

end subroutine res_rotation_precon
