!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!scale and shift interface points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine scale_shift_inter(x_inter)

  use interface_variables
  use fluid_variables,only:nsd

  real(8) x_inter(nsd,nn_inter)
  integer i

  do i=1,nn_inter
     x_inter(1:nsd,i)=x_inter(1:nsd,i)*scale_inter(1:nsd)
  end do
  do i=1,nn_inter
     x_inter(1:nsd,i)=x_inter(1:nsd,i)+shift_inter(1:nsd)
  end do

end subroutine scale_shift_inter

