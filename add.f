      subroutine add(newx, x, d, mn)

      implicit none
      include "global.h"

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
