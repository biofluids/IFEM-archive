      subroutine add(newx, x, d, mn)

      implicit none
      include "global.h"

      real* 8 x(mn,nnc), d(mn,nnc)
      real* 8 newx(mn,nnc)
      integer i,j,mn


      do i=1,nnc
         do j=1,mn
            newx(j,i) = x(j,i) + d(j,i)
         enddo
      enddo

      return
      end
