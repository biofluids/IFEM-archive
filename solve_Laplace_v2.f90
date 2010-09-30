!=============================================

!solve laplace eqn for the 2nd time

!============================================

subroutine solve_Laplace_v2(x_center,I_fluid_center,p_inter,w_inter,dg_inter,ne_outer,outer_ele,ne_inner,inner_ele,ien_center)

  use fluid_variables, only:nsd,ne
  use centermesh_variables, only:nn_center,ne_center,nen_center

  real(8) I_fluid_center(ne),x_center(nsd,ne)
  real(8) p_inter(nn_center),w_inter(nn_center),dg_inter(nn_center)
  integer ne_outer,ne_inner
  integer outer_ele(ne),inner_ele(ne)
  integer ien_center(nen_center,ne_center)

  integer i,j,icount,ie,inl,node
!...boundary condition
  I_fluid_center(:)=0.0
!  write(*,*)'inner_ele=',inner_ele(1:ne_inner)
  do ie=1,ne_outer
     I_fluid_center(outer_ele(ie))=0.0
  end do

  do ie=1,ne_inner
     I_fluid_center(inner_ele(ie))=1.0
  end do
!...initialize
  p_inter(:)=0.0
  w_inter(:)=0.0
!...calculate residual
  call block_Laplace(x_center,I_fluid_center,p_inter,w_inter,ien_center)
!...set id
  do ie=1,ne_outer
     p_inter(outer_ele(ie))=0.0
  end do

  do ie=1,ne_inner
     p_inter(inner_ele(ie))=0.0
  end do
!write(*,*)'p_inter=',p_inter(:)
  dg_inter(:)=0.0
  call gmres_Laplace_v2(x_center,I_fluid_center,w_inter,p_inter,dg_inter,ne_outer,outer_ele,ne_inner,inner_ele,ien_center)

  I_fluid_center(:)=I_fluid_center(:)+dg_inter(:)
!write(*,*)'I_fluid_center=',I_fluid_center(:)
end subroutine solve_Laplace_v2
  
