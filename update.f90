!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	updated.fcm                                                          c
!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine update(dinc, d, dg, mn)
  use fluid_variables
  implicit none

  integer,intent(in) :: mn
  real(8) :: d(mn,nn),dg(mn,nn)
  real(8),intent(out) :: dinc(mn,nn)

  dinc(1:mn,1:nn) = dg(1:mn,1:nn)

  d(1:mn,1:nn) = d(1:mn,1:nn) + dinc(1:mn,1:nn)

  return
end subroutine update
