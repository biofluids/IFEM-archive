c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  updated.fcm                                                          c
c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine update(d, dinc)

      implicit none
	  include "global.h"

	  real* 8 d(ndf,nnc), dinc(ndf,nnc)
	  integer idf, inl
	  
	  do inl = 1,nnc
		do idf = 1,ndf
		  d(idf,inl) = d(idf,inl) + dinc(idf,inl)
		end do
	  end do

	  return
	  end



