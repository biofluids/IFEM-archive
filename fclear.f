ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fclear(f,m)
      implicit none
      integer :: m
      real(8) :: f(m)

      integer :: i
      do i=1,m
         f(i) = 0.0
      enddo
      return
      end
