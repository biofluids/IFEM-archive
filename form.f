c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine formd(d)

      implicit none
	  include "global.h"

      real* 8  d(ndf,nnc)
	  integer idf,inn

      d(:,:) = 0.0
	  do inn = 1,nnc
		do idf = 1,ndf
		  d(idf,inn) = ic(idf)
		end do
	  end do
      
      return
      end
