!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	updated.fcm                                                          c
!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine update(dinc, d, dg, mn)
  use fluid_variables

  implicit none

  integer mn
  real* 8 d(mn,nn),dg(mn,nn)
  real* 8 dinc(mn,nn)
  !real* 8 hn(nn),hm(nn)
	
	
  !dinc(1:mn,1:nn) = 0.0d0

  !call equal(dg,dinc,mn*nn)
  dinc(1:mn,1:nn) = dg(1:mn,1:nn)

  !call gather(dinc, dg, mn, hn, hm)

  !do i=1,nn
  !   do j=1,mn 
  !      d(j,i) = d(j,i) + dinc(j,i)
  !   enddo
  !enddo

  d(1:mn,1:nn) = d(1:mn,1:nn) + dinc(1:mn,1:nn)

  return
end subroutine update
