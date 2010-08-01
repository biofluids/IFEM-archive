!
!  merge.f90
!  
!
!  Created by Xingshi Wang on 5/12/09.
!  Copyright 2009 __MyCompanyName__. All rights reserved.
!

       subroutine mergefinf(finf,nn_solid,mdata,n)
	   integer nn_solid
	   integer finf(nn_solid)
           integer finf1(nn_solid)
	   integer mdata(nn_solid)
	   integer n
	   integer tmparray(nn_solid)
	   integer i
	   integer j
	   integer bot
           finf1(:)=finf(:)
	   mdata(:)=0
	   bot = minval(finf1(1:nn_solid))
	   j=1
	   tmp=1
	   do while (tmp > 0) 
	          tmp = maxval(finf1(1:nn_solid))
			  do i=1,nn_solid
			     if (tmp .eq. finf1(i)) then
				     finf1(i)=0
				 endif
			   enddo
				 mdata(j)=tmp
				 j=j+1
	  enddo
	  n=j-2
			  
	   
			           
       
       return
       end
