subroutine r_nodalf
  use r_common
  use solid_variables, only: nn_solid
  implicit none

  integer :: i,ni


  if (numfn .gt. 0) then
     do i=1,numfn
        ni=i
        predrf(ni)=predrf(ni)+fnodo(i,1)
     enddo
  endif

  return
end subroutine r_nodalf
