!
!  getinf_el.f
!  
!
!  Created by Xingshi on 1/26/08.
!  Copyright 2008 __MyCompanyName__. All rights reserved.
!

       subroutine getinf_el_3d(finf, x, xna, nn, nsd, ne, nen, ien, maxconn)
	   implicit none
	   integer :: ninf, nn, nsd, ne, nen, maxconn ! fluid variable from main
	   integer ien(nen,ne) ! transfer matrix connectted local information to global
	   integer inf(maxconn) ! subroutine variable to store the index of element which contains the solid node
	   real(8) x(nsd), xna(nsd, nn) ! x - coordinates of solid node   xna - coordinates of fluid nodes in global view
           real(8) xn(nsd,nen)
           integer nnum
	   real(8) length_max(nsd)  ! the minimun rectangular to contian the element
	   real(8) length_min(nsd)
	   real(8) xna_p1(nsd)      ! coordinates of two neighbour points in the element
	   real(8) xna_p2(nsd)
	   real(8) xna_p3(nsd)
	   real(8) vector_0(nsd)
	   real(8) vector_1(nsd)
	   real(8) vector_2(nsd)
	   real(8) dot00
	   real(8) dot01
	   real(8) dot02
	   real(8) dot11
	   real(8) dot12
	   real(8) invDenom 
	   real(8) sign_1
	   real(8) sign_2   
	   real(8) vector_fe1(nsd)
	   real(8) vector_fe2(nsd)
	   real(8) vector_p(nsd)
	   real(8) sign_fe
	   integer finf  ! the final result we want to get
	   integer ie  ! loop variables
	   integer inl
	   integer isd  
	   integer i
       
       ninf = 0 
	   do ie = 1, ne
	      do inl = 1, nen

                     nnum=ien(inl,ie)
		         do isd = 1, nsd
			 xn(isd, inl) = xna(isd, nnum)
			 end do
	      end do
		 
		  do isd = 1, nsd
		  length_max(isd) = maxval(xn(isd,1:nen))
		  length_min(isd) = minval(xn(isd,1:nen))
                  end do
		
		  if (nsd==2) then
		     if (x(1)>length_min(1) .and. x(1)<length_max(1) .and. &
                         x(2)>length_min(2) .and. x(2)<length_max(2)) then
				 ninf = ninf +1
		if (ninf>maxconn) goto 2000
		         inf(ninf) = ie
			  end if
		   end if
		  
		    
		  if (nsd==3) then
		     if (x(1)>length_min(1) .and. x(1)<length_max(1) &
                         .and. x(2)>length_min(2) .and. x(2)<length_max(2) .and. &
				 x(3)>length_min(3) .and. x(3)<length_max(3)) then
				 ninf = ninf +1
		         inf(ninf) = ie
                        
			  end if
		   end if
           end do 
         !  write(*,*) 'ninf= ', ninf
                       
         !  write(*,*) 'inf= ' , inf(1)
         
	  if(ninf == 0) goto 2000
	   
	   if (nsd==2 .and. nen==3) then
	      do i= 1, ninf
				   xna_p1(1:nsd) = xna(1:nsd, ien(1,inf(i)))
				   xna_p2(1:nsd) = xna(1:nsd, ien(2,inf(i)))
				   xna_p3(1:nsd) = xna(1:nsd, ien(3,inf(i)))
				
			 vector_0(1:nsd) = xna_p3(1:nsd)-xna_p1(1:nsd)
			 vector_1(1:nsd) = xna_p2(1:nsd)-xna_p1(1:nsd)
			 vector_2(1:nsd) = x(1:nsd)-xna_p1(1:nsd)
			 
			 dot00=vector_0(1)*vector_0(1)+vector_0(2)*vector_0(2)
			 dot01=vector_0(1)*vector_1(1)+vector_0(2)*vector_1(2)
			 dot02=vector_0(1)*vector_2(1)+vector_0(2)*vector_2(2)
			 dot11=vector_1(1)*vector_1(1)+vector_1(2)*vector_1(2)
			 dot12=vector_1(1)*vector_2(1)+vector_1(2)*vector_2(2)
			 invDenom=1/(dot00*dot11-dot01*dot01)
			 sign_1=(dot11*dot02-dot01*dot12)*invDenom
			 sign_2=(dot00*dot12-dot01*dot02)*invDenom
			 
			 if(sign_1>0.0 .and. sign_2>0.0 .and. (sign_1+sign_2)<1.0) then
			 finf = inf(i)
			 goto 2000
			 end if
		  end do
           end if
          if (nsd==2 .and. nen==4) then
                 ! write(*,*) 'I am in the right loop' 
		  
                 ! write(*,*) 'get  x=' ,x(:)
                  do i= 1, ninf
                 ! write(*,*) 'xna= ',i, xna(1:nsd,ien(1:nen,inf(i)))
                  vector_fe1(1:nsd)=xna(1:nsd,ien(2,inf(i)))-xna(1:nsd,ien(1,inf(i)))
                  vector_fe2(1:nsd)=xna(1:nsd,ien(3,inf(i)))-xna(1:nsd,ien(1,inf(i)))
                  vector_p(1:nsd)=x(1:nsd)-xna(1:nsd,ien(1,inf(i)))
                  sign_1=vector_fe1(1)*vector_fe2(2)-vector_fe1(2)*vector_fe2(1)
                  sign_2=vector_fe1(1)*vector_p(2)-vector_fe1(2)*vector_p(1)
                  sign_fe=sign_2*sign_1
                 ! write(*,*) '1-2', sign_1, sign_2
                  if (sign_fe<0)  then
                  go to 2100
                  end if
                  !1-2
                  vector_fe1(1:nsd)=xna(1:nsd,ien(3,inf(i)))-xna(1:nsd,ien(2,inf(i)))
                  vector_fe2(1:nsd)=xna(1:nsd,ien(4,inf(i)))-xna(1:nsd,ien(2,inf(i)))
                  vector_p(1:nsd)=x(1:nsd)-xna(1:nsd,ien(2,inf(i)))
                  sign_1=vector_fe1(1)*vector_fe2(2)-vector_fe1(2)*vector_fe2(1)
                  sign_2=vector_fe1(1)*vector_p(2)-vector_fe1(2)*vector_p(1)
                  sign_fe=sign_2*sign_1
                 ! write(*,*) '2-3', sign_1, sign_2
                  if (sign_fe<0)  then
                  go to 2100
                  end if
                  !2-3
                  
                  vector_fe1(1:nsd)=xna(1:nsd,ien(4,inf(i)))-xna(1:nsd,ien(3,inf(i)))
                  vector_fe2(1:nsd)=xna(1:nsd,ien(1,inf(i)))-xna(1:nsd,ien(3,inf(i)))
                  vector_p(1:nsd)=x(1:nsd)-xna(1:nsd,ien(3,inf(i)))
                  sign_1=vector_fe1(1)*vector_fe2(2)-vector_fe1(2)*vector_fe2(1)
                  sign_2=vector_fe1(1)*vector_p(2)-vector_fe1(2)*vector_p(1)
                  sign_fe=sign_2*sign_1
                 ! write(*,*)'3-4',sign_1,sign_2
                  if (sign_fe<0)  then
                  go to 2100
                  end if
                  !3-4
                  
                  vector_fe1(1:nsd)=xna(1:nsd,ien(1,inf(i)))-xna(1:nsd,ien(4,inf(i)))
                  vector_fe2(1:nsd)=xna(1:nsd,ien(2,inf(i)))-xna(1:nsd,ien(4,inf(i)))
                  vector_p(1:nsd)=x(1:nsd)-xna(1:nsd,ien(4,inf(i)))
                  sign_1=vector_fe1(1)*vector_fe2(2)-vector_fe1(2)*vector_fe2(1)
                  sign_2=vector_fe1(1)*vector_p(2)-vector_fe1(2)*vector_p(1)
                  sign_fe=sign_2*sign_1
                 ! write(*,*) '4-1',sign_1,sign_2
                  if (sign_fe<0)  then
                  go to 2100
                  end if
                  !4-1
                  finf=inf(i)
                  go to 2000
2100              continue
                  end do
	  end if
		
	  if (nsd==3) then
!          write(*,*) 'i am here'
          call search_3d(finf, x, xna, nn, nsd, ne, nen, ien, inf, ninf,maxconn) ! call search_3d for 3d
          !cases
         end if	
		2000		continue
	!	write (*,*) 'Find the elemet in fluid domain as:',  finf
       return
       end
