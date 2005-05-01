!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c equal.f
!c let d=f, same dimensions
!ccccccccccccccccccccccccccccccccccccccccccccccc
subroutine equal(f,d,m)
      	implicit none
      
      	integer :: m
      	real(8) :: f(m),d(m)
		d = f
      	return
end
