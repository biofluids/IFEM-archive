!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!get boundary id for solving normal and curv!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_id_inter(ids,I_var)

  use fluid_variables, only:nn

  integer ids(nn)
  real(8) I_var(nn)

  integer i,j
  real(8) eps
  real(8) var1,var2

  ids(:) = 0
  eps = 1.0e-6

  do i=1,nn
     var1=abs(I_var(i))
     var2=abs(I_var(i)-1.0)
     if ((var1.le.eps).or.(var2.le.eps)) then
     ids(i) = 1
     end if
  end do

end subroutine get_id_inter

  
