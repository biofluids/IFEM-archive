subroutine getnorm(v1,v2,mn,norm)
  implicit none

  integer i,mn
  real*8 v1(mn),v2(mn),norm
  !real*8 mynorm,sdot
  !integer ierr

  norm = 0.0
  do i=1,mn
     norm =  norm + v1(i)*v2(i)
  enddo

  return
end
