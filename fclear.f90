!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fclear(f,m)
     	implicit none
      	integer :: m
      	real(8) :: f(m)

		f(:)=0.0
      	return
end
