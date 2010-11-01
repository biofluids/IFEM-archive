      subroutine add(newx, x, d, mn)
	use fluid_variables, only: nn
      implicit none

      real* 8 x(mn,nn), d(mn,nn)
      real* 8 newx(mn,nn)
      integer i,j,mn


      do i=1,nn
         do j=1,mn
            newx(j,i) = x(j,i) + d(j,i)
         enddo
      enddo

      return
      end
