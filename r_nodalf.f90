subroutine r_nodalf
  use r_common
  use solid_variables, only: nn_solid
  implicit none

  integer :: i,ni

  do i=1,numfn
     ni=(ndirfn(i)-1)*nn_solid+nodefn(i)
     predrf(ni)=predrf(ni)+fnodo(nodefn(i),ndirfn(i))
  enddo
  
  return
end subroutine r_nodalf
